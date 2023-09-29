# various tools
from dodo.parameters import AADICT_3_to_1, AADICT
from dodo.dodo_exceptions import dodoException
import matplotlib.pyplot as plt


# takes in lines from a pbd file and returns a dict.
def af2_lines_to_structure_dict(af2_pdb_lines):
    '''
    NOT USED ANYMORE!


    function to return a dictionary of the AF2 structure info.

    parameters
    ----------
    af2_pdb_lines : list
        the lines from the AF2 pdb file.

    Returns
    -------
    dict
        a dictionary of the structure info. Keys are:
        all_atom_coords : list
            a list of the coordinates of all atoms in the structure.
        ca_coords : list
            a list of the coordinates of all CA atoms in the structure.
        all_atoms_types : list
            a list of the atom types of all atoms in the structure.
        sequence : string
            the sequence of the structure in one letter code.
        ca_three_letter_seq : list
            a list of the three letter amino acid codes of all CA atoms in the structure.
        all_atom_three_letter_seq : list
            a list of the three letter amino acid codes of all atoms in the structure.
        coords_by_aa : dict
            For each amino acid, coords for every atom. Needed for the functionality to
            identify loops vs. IDRs
    '''
    # lists to hold info (and a str for sequence)
    all_atom_coords=[]
    ca_coords=[]
    all_atoms_types=[]
    sequence=''
    ca_three_letter_seq=[]
    all_atom_three_letter_seq=[]
    coords_by_aa_ind={}
    coords_by_aa_id=[]
    coord_by_aa_ind_with_name={}
    res_to_ind={}
    atom_coords_per_aa=[]
    start_resnum=1
    for line in af2_pdb_lines:
        if line != '':
            if line[:4]=='ATOM':
                resname=line[17:20]
                chain_id=line[21]                
                restype=line[13:16]
                restype_no_space=''
                for resval in restype:
                    if resval!=' ':
                        restype_no_space+=resval

                resnum=int(line[22:26])
                x=float(line[30:38])
                y=float(line[38:46])
                z=float(line[46:54])
                if resnum not in coords_by_aa_ind:
                    res_to_ind[resnum]=AADICT_3_to_1[resname]
                    coords_by_aa_ind[resnum]=[[x,y,z]]
                else:
                    cur_coords=coords_by_aa_ind[resnum]
                    cur_coords.append([x,y,z])
                    coords_by_aa_ind[resnum]=cur_coords
                if resnum not in coord_by_aa_ind_with_name:
                    coord_by_aa_ind_with_name[resnum]={restype_no_space:[x,y,z]}
                else:
                    cur_coord_by_aa_ind_with_name=coord_by_aa_ind_with_name[resnum]
                    cur_coord_by_aa_ind_with_name[restype_no_space]=[x,y,z]
                    coord_by_aa_ind_with_name[resnum]=cur_coord_by_aa_ind_with_name
                all_atom_coords.append([x,y,z])
                atom_coords_per_aa.append([x,y,z])
                all_atom_three_letter_seq.append(resname)
                all_atoms_types.append(restype_no_space)
                if restype_no_space=='CA':
                    ca_coords.append([x,y,z])
                    ca_three_letter_seq.append(resname)
                    sequence+=AADICT_3_to_1[resname]

    # popoulate the coords_by_aa_id dict
    atom_name_ind=0
    for aa in coords_by_aa_ind:
        cur_aa = res_to_ind[aa]
        cur_atom_coords=coords_by_aa_ind[aa]
        temp={}
        for a_ind in range(0, len(cur_atom_coords)):
            cur_atom_type=all_atoms_types[atom_name_ind]
            temp[cur_atom_type]=cur_atom_coords[a_ind]
            atom_name_ind+=1
        coords_by_aa_id.append({cur_aa:temp})


    # return lists
    return {'all_atom_coords':all_atom_coords, 'ca_coords':ca_coords, 
        'all_atoms_types':all_atoms_types, 'sequence':sequence, 
        'ca_three_letter_seq':ca_three_letter_seq, 
        'all_atom_three_letter_seq':all_atom_three_letter_seq,
        'coords_by_aa_ind':coords_by_aa_ind,
        'res_to_ind':res_to_ind,
        'coords_by_aa_id':coords_by_aa_id,
        'coord_by_aa_ind_with_name':coord_by_aa_ind_with_name,
        'final_coords_atoms_added':{}}





def plot_simple(all_coordinates, highlight=[], dims=[10,10]):
    '''
    simplified version of plot_structure. Doesn't require any coords
    info on fd vs. idr vs. loop.

    Parameters
    ----------
    all_coordinates : list
        a list of the coordinates of all atoms in the structure.
    highlight : list, optional
        a list of lists of coordinates to highlight. The default is [].
    dims : list, optional
        the dimensions of the plot. The default is [10,10].

    Returns
    -------
    None
        Just plots the graph and shows it.
    '''

    # Create a 3D plot
    fig = plt.figure(figsize=dims)
    ax = fig.add_subplot(111, projection='3d')
    if highlight != []:
        for coord in highlight:
            all_coordinates.remove(coord)
    x_c, y_c, z_c = zip(*all_coordinates)
    ax.scatter(x_c, y_c, z_c, s=70, depthshade=True)
    if highlight != []:
        hi_xc, hi_yc, hi_zc=zip(*highlight)
        for hi in range(0, len(highlight)):
            curx=hi_xc[hi]
            cury=hi_yc[hi]
            curz=hi_zc[hi]
            ax.scatter(curx, cury, curz, s=70, depthshade=True, color='red')
    
    # Set make sure axes all same size so structure looks good. 
    all_coords=[]
    all_coords.extend(x_c)
    all_coords.extend(y_c)
    all_coords.extend(z_c)
    ax.set_ylim(min(all_coords),max(all_coords))
    ax.set_xlim(min(all_coords),max(all_coords))
    ax.set_zlim(min(all_coords),max(all_coords))
    
    # get rid of annoying stuffs
    all_ax=[ax.xaxis, ax.yaxis, ax.zaxis]
    for curax in all_ax:
        curax.set_pane_color((1.0, 1.0, 1.0, 0.0))
        curax.set_ticks([])
        curax._axinfo["grid"]['color'] =  (1,1,1,0)
        curax._axinfo['axisline']['linewidth']=0
        curax._axinfo['grid']['linewidth']=0
        curax._axinfo['juggled']=(0,0,0)

    plt.show()


def plot_structure(all_coordinates, fd_idr_loop_dict, dims=[10,10], save_path='', 
    fd_color='blue', IDR_color='orange', loop_color='red'):
    '''
    Quick way to view a structure using matplotlib.

    Parameters
    ----------
    all_coordinates : list
        a list of the coordinates of all atoms in the structure.
    fd_idr_loop_dict : dict
        a dictionary of the structure info. keys have 'idr'
        in them for idrs, 'folded' in them for FDs, and
        'loop' in them for loops. Values are lists 
        that contain the coordinates in the sequence for each
        different 'feature' type.
    dims : list, optional
        the dimensions of the plot. The default is [10,10].
    save_path : string, optional
        if you want to save the plot, specify a path here.
        The default is ''. If it is just not specified, just shows plot.
    fd_color : string, optional
        the color to plot the folded domains. The default is 'blue'.
    IDR_color : string, optional
        the color to plot the IDRs. The default is 'orange'.
    loop_color : string, optional
        the color to plot the loops. The default is 'red'.


    Returns
    -------
    None
        Shows plot of the structure.
    '''
    # Create a 3D plot
    fig = plt.figure(figsize=dims)
    ax = fig.add_subplot(111, projection='3d')

    # break up coords into each 'part'
    fd_coords=[]
    idr_coords=[]
    loop_coords=[]
    for fd_idr_loop in fd_idr_loop_dict:
        if 'idr' in fd_idr_loop:
            cur_coords=fd_idr_loop_dict[fd_idr_loop]
            idr_coords.append(all_coordinates[cur_coords[0]:cur_coords[1]])
        elif 'folded' in fd_idr_loop:
            cur_coords=fd_idr_loop_dict[fd_idr_loop]
            fd_coords.extend(all_coordinates[cur_coords[0]:cur_coords[1]])
        elif 'loop' in fd_idr_loop:
            fd_loop_coords=fd_idr_loop_dict[fd_idr_loop]
            cur_fd_coords=fd_loop_coords[0]
            cur_loop_in_fd_coords=fd_loop_coords[1]
            all_fd_coords=all_coordinates[cur_fd_coords[0]:cur_fd_coords[1]]
            for loop in cur_loop_in_fd_coords:

                cur_loop_coords=all_coordinates[loop[0]:loop[1]]
                for coord in cur_loop_coords:
                    all_fd_coords.remove(coord)
                loop_coords.append(cur_loop_coords)
            fd_coords.extend(all_fd_coords)
        else:
            raise dodoException('Error - coordinate input did not seem to be an IDR, FD, or loop.')

    # extract coords for each.
    if fd_coords != []:
        fd_x_coords, fd_y_coords, fd_z_coords = zip(*fd_coords)
        if fd_color != 'blue':
            ax.scatter(fd_x_coords, fd_y_coords, fd_z_coords, s=70, depthshade=True, color=fd_color,)
        else:
            ax.scatter(fd_x_coords, fd_y_coords, fd_z_coords, s=70, depthshade=True)


    if idr_coords != []:
        if len(idr_coords) > 0:
            for cur_coord in idr_coords:
                idr_x_coords, idr_y_coords, idr_z_coords = zip(*cur_coord)
                # plt the IDR
                ax.plot(idr_x_coords, idr_y_coords, idr_z_coords, color=IDR_color, linewidth=2) 
        else:   
            idr_x_coords, idr_y_coords, idr_z_coords = zip(*idr_coords)
            ax.plot(idr_x_coords, idr_y_coords, idr_z_coords, color=IDR_color, linewidth=2) 

    if loop_coords != []:
        if len(loop_coords) > 1:
            for cur_coord in loop_coords:
                loop_x_coords, loop_y_coords, loop_z_coords = zip(*cur_coord)  
                ax.plot(loop_x_coords, loop_y_coords, loop_z_coords, color=loop_color, linewidth=2)              
        else:
            loop_x_coords, loop_y_coords, loop_z_coords = zip(*loop_coords[0])
            # plot the loops
            ax.plot(loop_x_coords, loop_y_coords, loop_z_coords, color=loop_color, linewidth=2)



    # Set make sure axes all same size so structure looks good. 
    all_coords_x, all_coords_y, all_coords_z = zip(*all_coordinates)
    all_coords=[]
    all_coords.extend(all_coords_x)
    all_coords.extend(all_coords_y)
    all_coords.extend(all_coords_z)
    ax.set_ylim(min(all_coords),max(all_coords))
    ax.set_xlim(min(all_coords),max(all_coords))
    ax.set_zlim(min(all_coords),max(all_coords))
    
    # get rid of annoying stuffs
    all_ax=[ax.xaxis, ax.yaxis, ax.zaxis]
    for curax in all_ax:
        curax.set_pane_color((1.0, 1.0, 1.0, 0.0))
        curax.set_ticks([])
        curax._axinfo["grid"]['color'] =  (1,1,1,0)
        curax._axinfo['axisline']['linewidth']=0
        curax._axinfo['grid']['linewidth']=0
        curax._axinfo['juggled']=(0,0,0)

    # save or show
    if save_path=='':
        # Show the plot
        plt.show()
    else:
        plt.savefig(save_path)

