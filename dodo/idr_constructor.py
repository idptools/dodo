# imports from 'the web' as the kids say
import os 
import random
import math
import numpy as np
from copy import deepcopy, copy
from sparrow import Protein as pr

# dodo imports
from dodo.pdb_tools import PDBParser
# for the ... problems.
from dodo.dodo_exceptions import dodoException
from dodo import parameters
from dodo.all_atoms import add_necessary_C_N

def get_res_dist(xyz1, xyz2):
    '''
    get distance between two coordinates in 3D space

    Parameters
    ----------
    xyz1 : list
        list of x,y,z coordinates for first point
    xyz2 : list
        list of x,y,z coordinates for second point

    Returns
    -------
    float
        distance between two points
    '''
    xdiff = xyz1[0]-xyz2[0]
    ydiff = xyz1[1]-xyz2[1]
    zdiff = xyz1[2]-xyz2[2]
    return math.sqrt((xdiff*xdiff)+(ydiff*ydiff)+(zdiff*zdiff))


def is_point_within_sphere_constraints(point, center, radius, off_by=None, verbose=False):
    """
    Determine whether a 3D coordinate is within a sphere.
    What we can do with this is basically have an objective end point
    in our folded domain we want to connect the IDR to and shrink that 
    radius over time so we eventually get ~3.8Å from the final position.

    Parameters
    ----------
    point : (tuple or list): 
        The 3D coordinate (x, y, z) to check.
    center : (tuple or list)
        The center of the sphere (x_center, y_center, z_center).
    radius : (float)
        The radius of the sphere.

    Returns
    ----------
    bool : 
        True if the point is inside the sphere, False otherwise.
    """
    if len(point) != 3 or len(center) != 3:
        raise ValueError("Both point and center must be 3D coordinates.")

    # Calculate the Euclidean distance between the point and the sphere's center
    distance = math.sqrt((point[0] - center[0])**2 + (point[1] - center[1])**2 + (point[2] - center[2])**2)
    
    # Check if the distance is less than or equal to the sphere's radius
    if off_by != None:
        if verbose:
            print(f'objective min: {radius-off_by}, objective_max: {radius} actual: {distance}')
        if distance <= radius and distance > radius-off_by:
            return True
    else:
        return distance <= radius


# should make this more stream lined but honestly don't care so :shrug emoji SRY JEFF
def remove_IDRs_loops(PDBParserObj):
    '''
    Removes idrs and loops so we can move stuff around and then
    rebuild. I did it in a stupid way before like a champ.
    
    Parameters
    -----------
    PDBParserObj
        PDBParser object

    Returns
    --------
    PDBParserObj
    '''
    seq=PDBParserObj.sequence
    nuke_these_coords = []
    # get loop vals
    for loop in PDBParserObj.FD_loop_coords:
        cur_loop=PDBParserObj.FD_loop_coords[loop]
        for val in cur_loop[1]:
            loop_coords=[aa for aa in range(val[0], val[1])]
            nuke_these_coords.extend(loop_coords)
    for idr in PDBParserObj.IDR_coords:
        cur_idr = PDBParserObj.IDR_coords[idr]
        # make sure we get the final value...
        if cur_idr[-1]==len(seq)-1:
            nuke_these_coords.extend([aa for aa in range(cur_idr[0], cur_idr[1]+1)])    
        else:
            nuke_these_coords.extend([aa for aa in range(cur_idr[0], cur_idr[1])])
    # NUKE.
    for aa in nuke_these_coords:
        PDBParserObj.all_atom_coords_by_index[aa] = {} 
    return PDBParserObj

def overwrite_beta_for_vis(PDBParserObj, FD=100, IDR=0):
    '''
    Overwrites the beta values in a PDBParserObj
    Can use for visualization purpose in VMD.
    
    Parameters
    -----------
    PDBParserObj
        PDBParser object
    FD : int
        FD number to overwrite with
        default: 0
    IDR : int
        IDR number to overwrite with
        default: 100

    Returns
    --------
    PDBParserObj
    '''
    seq=PDBParserObj.sequence
    IDR_coords_to_overwrite = []
    # get loop vals
    for loop in PDBParserObj.FD_loop_coords:
        cur_loop=PDBParserObj.FD_loop_coords[loop]
        for val in cur_loop[1]:
            loop_coords=[aa for aa in range(val[0], val[1])]
            IDR_coords_to_overwrite.extend(loop_coords)
    for idr in PDBParserObj.IDR_coords:
        cur_idr = PDBParserObj.IDR_coords[idr]
        # make sure we get the final value...
        if cur_idr[-1]==len(seq)-1:
            IDR_coords_to_overwrite.extend([aa for aa in range(cur_idr[0], cur_idr[1]+1)])    
        else:
            IDR_coords_to_overwrite.extend([aa for aa in range(cur_idr[0], cur_idr[1])])
    # overwrite.
    for res in PDBParserObj.beta_vals_by_index:
        if res in IDR_coords_to_overwrite:
            PDBParserObj.beta_vals_by_index[res]=IDR
        else:
            PDBParserObj.beta_vals_by_index[res]=FD
    return PDBParserObj


def get_current_alpha_carbons(PDBParserObj, return_list=True):
    '''
    Quick function to get the current alpha carbons from the PDBParserObj.

    Parameters
    ----------
    PDBParserObj : PDBParser
        PDBParser object
    return_list : bool
        if True, returns list of coordinates. If False, returns dict of coordinates.
        default: True
    '''
    # get all atom
    all_atom_coords=PDBParserObj.all_atom_coords_by_index
    # get dict
    fin_dict={}
    for aa in all_atom_coords:
        if 'CA' in all_atom_coords[aa].keys():
            fin_dict[aa]=all_atom_coords[aa]['CA']
    if return_list:
        return list(fin_dict.values())
    else:
        return fin_dict


def get_region_coords(PDBParserObj, specific_regions=[[]], return_list=True):
    '''
    Handy function to get the coordinates of a specific region.
    Can specify multiple regions.

    Parameters
    ----------
    PDBParserObj : PDBParser
        PDBParser object
    specific_regions : list of lists of regions. If empty, returns all atom coords.
        default: [[]]
    return_list : bool
        if True, returns list of coordinates. If False, returns dict of coordinates.
        default: True

    Returns
    -------
    list or dict
        list or dict of coordinates
    '''
    if type(specific_regions[0])!= list:
        raise dodoException('if you specify regions, be sure to input a list of lists!')
    # get all atom
    all_atom_coords=PDBParserObj.all_atom_coords_by_index
    # get regions.
    fin_vals=[]
    if specific_regions==[[]]:
        coords=list(all_atom_coords.keys())
    else:
        coords=[]
        for reg in specific_regions:
            coords.extend([aa for aa in range(reg[0], reg[1]+1)])
    # get dict
    fin_dict={}
    for aa in coords:
        if aa not in fin_dict:
            fin_dict[aa]={}
        for at in all_atom_coords[aa]:
            curfin=fin_dict[aa]
            curfin[at]=all_atom_coords[aa][at]
            fin_dict[aa]=curfin
    if return_list:
        return list(fin_dict.values())
    else:
        return fin_dict



def is_min_distance_satisfied(PDBParserObj, coordinate_of_interest_ind, 
    coordinate_of_interest_coordinate, coordinate_of_interest_atom,
    average_bonds = parameters.average_bond_distances_AF2,
    CA_clash_dist=parameters.CA_clash_dist, bond_buffer=parameters.bond_buffer,
    min_dist_adjacent=parameters.min_dist_adjacent, verbose=False):
    '''
    Function to see if a specific coordinate is within a distance of all
    cooridnates in a list of coordinates.

    Parameters
    ----------
    PDBParserObj : PDBParser
        PDBParser object
    coordinate_of_interest_ind : int
        index of coordinate of interest that we are checking for clashes against.
    coordinate_of_interest_coordinate : list
        list of x,y,z coordinates for coordinate of interest
    coordinate_of_interest_atom : string
        atom name of coordinate of interest 
        ex. 'CA'
    average_bonds : dict
        dict of average bond distances. Default is those specified in parameters.py
    CA_clash_dist : float
        minimum distance between CA atoms. Default is those specified in parameters.py
    bond_buffer : float
        buffer to add to average bond distances. Default is those specified in parameters.py
        This lets us have things a little closer than the bond length but not by much.
    min_dist_adjacent : float
        minimum distance between adjacent atoms. Default is those specified in parameters.py
    verbose : bool
        if True, prints out why min distance was not satisfied. Default is False.

    Returns
    -------
    bool
        True if min distance is satisfied, False otherwise.
    '''
    # Get all atom coords
    all_atom_coords = PDBParserObj.all_atom_coords_by_index
    # get rid of empty dicts
    all_atom_coords = {aa:atom for aa,atom in all_atom_coords.items() if atom != {}}
    # get immediately adjacent atom indices to coord of interest
    if coordinate_of_interest_ind !=0:
        surround_coordinates_and_atoms=[coordinate_of_interest_ind-1,coordinate_of_interest_ind+1]
    else:
        surround_coordinates_and_atoms=[coordinate_of_interest_ind+1]
    
    # get final coords to check
    final_coords_not_adjacent=[]
    adjacent_coords={}
    for aa in all_atom_coords:
        if aa not in surround_coordinates_and_atoms:
            final_coords_not_adjacent.extend(list(all_atom_coords[aa].values()))
        else:
            adjacent_coords[aa]=all_atom_coords[aa]
    
    # use the CA clash dist for min distances of non-adjacent
    # first make into arrays
    final_coords_not_adjacent=np.array(final_coords_not_adjacent)
    coord_of_interest_array = np.array(coordinate_of_interest_coordinate)

    # Calculate the Euclidean distances between the coordinate of interest and all points
    distances_non_adjacent = np.linalg.norm(final_coords_not_adjacent - coord_of_interest_array, axis=1)

    # if not min CA distances, return False
    if np.all(distances_non_adjacent >= CA_clash_dist)!=True:
        if verbose:
            print('CA_fail')
        return False

    # If we pass the CA checks, now check adjacent atoms.
    # this checks all atoms of adjacent atoms. These distances 
    # that are specified for the clash are in parameters.py for 
    #each possible pair of atoms
    adjacent_atom_indices=list(adjacent_coords.keys())
    if len(adjacent_atom_indices)==2:
        prev_adj_atoms=adjacent_coords[adjacent_atom_indices[0]]
        for at in prev_adj_atoms:
            at_pair = f'{at}_CA'
            if at_pair in average_bonds:
                # get approx min dist. 
                min_dist=average_bonds[at_pair]-0.1
                if get_res_dist(prev_adj_atoms[at], coordinate_of_interest_coordinate)< min_dist:
                    if verbose:
                        print('pre_adj_atom fail')
                    return False
            else:
                if get_res_dist(prev_adj_atoms[at], coordinate_of_interest_coordinate)< min_dist_adjacent:
                    if verbose:
                        print('pre_adj_atom fail')
                    return False
    # check next atom
    if len(adjacent_coords)==1:
        next_adj_atoms=adjacent_coords[adjacent_atom_indices[0]]
        for at in next_adj_atoms:
            at_pair = f'CA_{at}'
            if at_pair in average_bonds:
                # get approx min dist. 
                min_dist=average_bonds[at_pair]-0.1
                if get_res_dist(next_adj_atoms[at], coordinate_of_interest_coordinate)< min_dist:
                    if verbose:
                        print('post_adj_atom fail')                
                    return False
            else:
                if get_res_dist(next_adj_atoms[at], coordinate_of_interest_coordinate)< min_dist_adjacent:
                    if verbose:
                        print('pre_adj_atom fail')                
                    return False  
    # if get down to here, return True
    return True  


def is_min_distance_satisfied_idr_from_seq(coords_list, coordinate_of_interest, 
    min_distance=parameters.CA_clash_dist):
    '''
    Function to see if a specific coordinate is within a distance of all
    cooridnates in a list of coordinates.

    Parameters
    ----------
    coords_list : list
        list of coordinates
    coordinate_of_interest : tuple
        tuple of x,y,z coordinates
    min_distance : float
        minimum distance to check

    Returns
    -------
    bool
        True if min distance is satisfied, False otherwise.
    '''

    # Convert the input lists to NumPy arrays
    coords_array = np.array(coords_list)
    coord_of_interest_array = np.array(coordinate_of_interest)

    # Calculate the Euclidean distances between the coordinate of interest and all points
    distances = np.linalg.norm(coords_array - coord_of_interest_array, axis=1)

    # Check if the minimum distance is greater than or equal to the specified min_distance
    return np.all(distances >= min_distance)


def is_min_dist_satisfied_FDs(PDBParserObj, moved_fd_coordinates_indices, 
    check_these_coords=[[]], CA_clash_dist=parameters.CA_clash_dist):
    '''
    function to check for clashes between fds after moiving. Difference here
    is that we are just checking the CA coords. This should be sufficient for the FDs.

    Parameters
    ----------
    PDBParserObj : PDBParser
        PDBParser object
    moved_fd_coordinates_indices : list
        list of indices for the FD that was moved
    check_these_coords : list of lists
        list of indices for the FDs that were already moved. Don't check FDs we haven't moved yet!
        default: [[]]
        This lets us not check the coordinates of FDs before we move them.
    CA_clash_dist : float
        minimum distance between CA atoms. Default is those specified in parameters.py

    Returns
    -------
    bool
        True if min distance is satisfied, False otherwise.
    '''
    # get coords for moved fd. 
    moved_fd = get_region_coords(PDBParserObj, specific_regions=[moved_fd_coordinates_indices], return_list=True)

    # get list of coords not in fd. 
    all_coords = get_region_coords(PDBParserObj, return_list=False, specific_regions=check_these_coords)

    # make coords from that into list, provided they aren't in the all_
    all_other_coords=[]
    all_moved_fd_coord_indices=[aa for aa in range(moved_fd_coordinates_indices[0], moved_fd_coordinates_indices[1])]
    for aa in all_coords:
        if aa not in all_moved_fd_coord_indices:
            for at in all_coords[aa]:
                all_other_coords.append(all_coords[aa][at])
    # now check for clashes
    for aa in moved_fd:
        for at in aa:
            coord=aa[at]
            # use the CA clash dist for min distances of non-adjacent
            # first make into arrays
            all_other_coords=np.array(all_other_coords)
            coord_of_interest_array = np.array(coord)
            # Calculate the Euclidean distances between the coordinate of interest and all points
            distances = np.linalg.norm(all_other_coords - coord_of_interest_array, axis=1)
            # if not min CA distances, return False
            if np.all(distances >= CA_clash_dist)!=True:
                return False
    # otherwise return True
    return True


def random_coordinate_on_sphere_surface(center_coordinate, radius):
    """
    Generate a random 3D coordinate on the surface of a sphere with a specified radius.

    Parameters
    ----------
    center_coordinate : tuple
        The center of the sphere (x_center, y_center, z_center).
    radius : float
        The radius of the sphere.

    Returns
    --------
    tuple
        A random 3D coordinate on the sphere's surface.
    """
    if len(center_coordinate) != 3:
        raise ValueError("Center coordinate must be a 3D coordinate (x_center, y_center, z_center).")

    if radius <= 0:
        raise ValueError("Radius must be greater than 0.")

    # Generate random spherical coordinates
    theta = random.uniform(0, 2 * math.pi)  # Azimuthal angle
    phi = random.uniform(0, math.pi)  # Polar angle

    # Convert spherical coordinates to Cartesian coordinates
    x = center_coordinate[0] + radius * math.sin(phi) * math.cos(theta)
    y = center_coordinate[1] + radius * math.sin(phi) * math.sin(theta)
    z = center_coordinate[2] + radius * math.cos(phi)

    return (x, y, z)

def create_point_on_line(start, end, distance):
    '''
    Function to create a point on a line between two points.
    
    Parameters
    ----------
    start : tuple
        tuple of x,y,z coordinates for start point
    end : tuple
        tuple of x,y,z coordinates for end point
    distance : float
        distance to move from start point

    Returns
    -------
    tuple
        tuple of x,y,z coordinates for new point
    '''
    # Convert the input coordinates to numpy arrays for vector operations
    start = np.array(start)
    end = np.array(end)

    # Calculate the direction vector from start to end
    direction_vector = end - start

    # Normalize the direction vector to have a unit length
    direction_unit_vector = direction_vector / np.linalg.norm(direction_vector)

    # Calculate the new point by adding the scaled direction vector to the start point
    new_point = start + (distance * direction_unit_vector)

    return new_point.tolist()


def predict_e2e(sequence): 
    '''
    function to predict end to end distance using SPARROW. 
    Just a wrapper because I'm too lazy to type the thing every time
    and I always forget exactly how to do it because I'm a few Crayons™ 
    short of a full box.

    Parameters
    ----------
    sequence : string
        the amino acid sequence as a string. 

    Returns
    --------
    float
        the predicted end to end distance in Å.
    '''
    return pr(sequence).predictor.end_to_end_distance()


def calculate_center_coordinate(coordinates):
    """
    Calculate the ~center coordinate of a list of 3D coordinates.
    This isn't COM. Rather, we just need an approximate center of the FDs
    so we can draw a line between those centers and translate the FDs relative to
    each other. 

    Parameters
    -----------
        coordinates : (list of tuples) 
            List of 3D coordinates [(x1, y1, z1), (x2, y2, z2), ...].

    Returns
    -------
    tuple
        The center coordinate (x_center, y_center, z_center).
    """
    if not coordinates:
        raise ValueError("List of coordinates is empty.")

    num_coordinates = len(coordinates)
    x_sum, y_sum, z_sum = 0, 0, 0

    for coord in coordinates:
        if len(coord) != 3:
            raise ValueError("All coordinates must be 3D, i.e., have three components (x, y, z).")
        x, y, z = coord
        x_sum += x
        y_sum += y
        z_sum += z

    x_center = x_sum / num_coordinates
    y_center = y_sum / num_coordinates
    z_center = z_sum / num_coordinates

    return (x_center, y_center, z_center) 

def translate_coordinates(coordinates, specific_coordinate, target_coordinate):
    """
    Translate all coordinates in a list such that a specific coordinate is 
    located at a target coordinate. This function lets us take the folded 
    domains and move them while holding the actual integrity (defined
    as the coordinates relative to each other in the domain) constant. 
    Specifically, we can move the FDs to be closer or further away, 
    letting us build IDRs that are more expanded or compact between them.

    Parameters
    -----------
    coordinates : (list of tuples)
        List of 3D coordinates [(x1, y1, z1), (x2, y2, z2), ...].
    specific_coordinate : (tuple) 
        The specific 3D coordinate to be moved.
    target_coordinate : (tuple) 
        The 3D coordinate where the specific coordinate should be placed.

    Returns
    -----------
    list of tuples : 
        List of translated 3D coordinates.
    """
    if len(specific_coordinate) != 3 or len(target_coordinate) != 3:
        raise ValueError("Both specific_coordinate and target_coordinate must be 3D coordinates.")

    # Calculate the translation vector
    translation_vector = np.array(target_coordinate) - np.array(specific_coordinate)

    # Translate all coordinates in the list
    translated_coordinates = [(x + translation_vector[0], y + translation_vector[1], z + translation_vector[2])
                              for x, y, z in coordinates]

    return translated_coordinates


def generate_coordinate_on_line(coord1, coord2, distance, from_coord1=True):
    """
    Generate a 3D coordinate along a straight line, given two coordinates 
    that is a specific distance from the first coordinate.
    This lets us translate the FDs along a specific line. Specifically, we can
    take the first FD and the second FD, get their centers (approximately), 
    draw a line between them, and then translate them in relative 3D space along
    a specific line using the translate_coordinates function.

    Parameters
    ----------
    coord1 : tuple
        The first 3D coordinate (x1, y1, z1).
    coord2 : tuple
        The second 3D coordinate (x2, y2, z2).
    distance : float
        The specified distance between the coordinates.
    from_coord1 : bool
        whether the returned coordinate should be the specified distance
        from coord1 or from coord2.
        Default is coord1. If False, distance is from coord2
    Returns
    -------
    list of tuples : 
        List of 3D coordinates forming a straight line.
    """
    # check those coords.
    if len(coord1) != 3 or len(coord2) != 3:
        raise ValueError("Both coordinates must be 3D, i.e., have three components (x, y, z).")

    # Convert coordinates to numpy arrays for vector operations
    coord1 = np.array(coord1)
    coord2 = np.array(coord2)

    # Calculate the direction vector between coord1 and coord2
    if from_coord1:
        direction_vector = coord2 - coord1
    else:
        direction_vector = coord1 - coord2

    # Calculate the unit vector in the direction of the line
    unit_vector = direction_vector / np.linalg.norm(direction_vector)

    # Generate coordinate along the line at specified interval
    if from_coord1:
        return coord1+(distance*unit_vector)
    else:
        return coord2+(distance*unit_vector)


def round_coordinates(PDBParserObj, decimals=3):
    '''
    quick function to do final rounding after all coords are generated.
    Going to avoid rounding before final generation because translation of
    FDs is somewhat iterative, so rounding iteratively could compound the
    error of multi FD proteins. 

    Parameters
    ----------
    PDPParserObj
        PDBParser object

    decimals : int
        number of decimals to round to
        default: 3

    Returns
    -------
    PDBParserObj
        returns PDBParserObj with all coords rounded. 
    '''
    # get all coordinates
    all_coords=PDBParserObj.all_atom_coords_by_index

    # iterate through coordinates
    for aa in all_coords:
        for at in all_coords[aa]:
            # get_coords
            cur_coords=np.array(all_coords[aa][at])
            
            x=round(cur_coords[0], decimals)
            y=round(cur_coords[1], decimals)
            z=round(cur_coords[2], decimals)
                # update dict
            PDBParserObj.all_atom_coords_by_index[aa][at]=(x,y,z)
    return PDBParserObj

def random_walk_step(initial_coords, distance):
    '''
    Function to talk a random walk step.
    
    Parameters
    ----------
    initial_coords : tuple
        tuple of x,y,z coordinates for initial point
    distance : float
        distance to move from initial point

    Returns
    -------
    tuple
        tuple of x,y,z coordinates for new point
    '''
    # maths stuffs
    theta = random.uniform(0, 2 * math.pi)
    phi = random.uniform(0, math.pi)

    # Calculate the Cartesian coordinates
    x = initial_coords[0] + distance * math.sin(phi) * math.cos(theta)
    y = initial_coords[1] + distance * math.sin(phi) * math.sin(theta)
    z = initial_coords[2] + distance * math.cos(phi)

    return (x, y, z)


def randomish_step(start_coordinate, final_coordinate, step_outside_final_coord,
    step_from_start_coord_dist=3.8):
    '''
    By taking a random walk step from the final coordinate, we can
    basically move in that direction. 

    Parameters
    ----------
    start_coordinate : tuple
        tuple of x,y,z coordinates for start point
    final_coordinate : tuple
        tuple of x,y,z coordinates for end point
    step_outside_final_coord : float
        distance from outside final coordinate to move to
        higher numbers = less likely to be right at the cooridnate. 
    step_from_start_coord_dist : float
        the distance from the start coordiante to ultiimately traverse

    Returns
    -------
    tuple
        tuple of x,y,z coordinates for new point\
    '''
    # randomly choose another direction.
    random_direction=random_walk_step(final_coordinate, step_outside_final_coord)
    return generate_coordinate_on_line(start_coordinate, random_direction, 
        step_from_start_coord_dist, from_coord1=True)


def get_seq_from_all_atom_coords(PDBParserObj):
    '''
    Function to get the sequence from the PDBParserObj
    using the all atom coordinates. Because the all atom
    coordinates are updated as we build the various domains,
    this lets us do a final check that the final sequence we build
    is the expected sequence. 

    Parameters
    ----------
    PDBParserObj : PDBParser
        PDBParser object

    Returns
    -------
    string
        string of amino acids
    '''
    # get all atom coords
    all_atom_coords=PDBParserObj.all_atom_coords_by_index
    # get the amino acids corresponding to each index.
    index_to_aa=PDBParserObj.index_to_3aa
    # empty string to hold sequence
    seq=''
    # iterate through by each aa. 
    for aa_ind in all_atom_coords:
        amino_acid=index_to_aa[aa_ind]
        seq+=parameters.AADICT_3_to_1[amino_acid]
    return seq

def center_fd_at_0(PDBParserObj):
    '''
    function to take the first fd and center it at 000.
    '''
    # get atoms
    atoms=PDBParserObj.all_atom_coords_by_index
    # get those FDs
    get_fds=PDBParserObj.FD_coords
    # get loops
    fd_loops=PDBParserObj.FD_loop_coords
    for loop_num in fd_loops:
        get_fds[loop_num]=(fd_loops[loop_num][0])
    all_fds={}
    for fd in get_fds:
        fd_num=int(fd.split('_')[-1])
        all_fds[fd_num]=get_fds[fd]
    # sort
    fds_in_order=sorted(all_fds)
    # get coords of first fd
    domain1_coords_ind = all_fds[fds_in_order[0]]
    # get the ref point of the final CA of the first fd, we move this to 000
    ref_point1 = np.array(atoms[domain1_coords_ind[-1]-1]['CA'])   
    # Get actual coords for fd regions.
    domain1_coords=get_region_coords(PDBParserObj, specific_regions=[domain1_coords_ind], return_list=False)            
    
    # for repopulating PDBParser after confirming everything is good..
    aa_and_at_dom1=[]
    domain1_coords_list=[]
    for aa in domain1_coords:
        for at in domain1_coords[aa]:
            aa_and_at_dom1.append([aa, at])
            domain1_coords_list.append(domain1_coords[aa][at])

    # now have list of coords, need to translate.
    # translate coords of fd to 000
    translated_coords=translate_coordinates(domain1_coords_list, ref_point1, (0,0,0))
    #IDK why I did it this way but ... yup
    # track numbers of aa and atom for repopulating dict
    raw_num_ind=0
    for coord in translated_coords:
        # get aa num and atom name to repopulate thigns
        cur_pair = aa_and_at_dom1[raw_num_ind]
        cur_ats = PDBParserObj.all_atom_coords_by_index[cur_pair[0]]
        cur_ats[cur_pair[1]]=coord
        # update all atom coords dict in pdbparserobj
        PDBParserObj.all_atom_coords_by_index[cur_pair[0]]=cur_ats
        # update the index so we track atoms and aa numbers
        raw_num_ind+=1
    return PDBParserObj



def place_FDs(PDBParserObj, mode, linear_placement=False,
    verbose=True, total_attempts=30):
    '''
    Function to place the FDs. Basically iterates through all the FDs and 
    FDs with loops and places them at random locations at specific distances
    that are dependent on mode. After this, we can rebuild the IDRs. 

    Parameters
    ----------
    PDBParserObj : PDBParser
        PDBParser object
    mode : string
        mode to use for placing FDs. Can be super_compact, compact, normal, 
        expanded, super_expanded, max_expansion, or predicted. 
        Predicted predicts the end-to-end distance of each IDR and places the FDs
        such that they will be that distance away from each other. 
    linear_placement : bool
        whether to place the folded domains across a linear axis. 
        Default : False
    verbose : bool
        If True, prints out what is happening. Default is True.
    total_attempts : int
        total number of attempts to try to place FDs. Default is 30.

    Returns
    -------
    PDBParserObj
        PDBParser object with FDs placed.
    '''

    # check the mode
    run_prediction=False
    if mode =='super_compact':
        length_multiplier = 0.18
    elif mode == 'compact':
        length_multiplier=0.55
    elif mode == 'normal':
        length_multiplier=0.8
    elif mode == 'expanded':
        length_multiplier = 1.05
    elif mode == 'super_expanded':
        length_multiplier = 1.3
    elif mode == 'max_expansion':
        length_multiplier=1.65
    elif mode == 'predicted':
        run_prediction=True
    else:
        raise dodoException('Invalid mode specified. Please specify super_compact, compact, normal, expanded, super_expanded, or max_expansion, or predicted.')


    # now we need to move the FDs relative to each other. Originally I tried to be clever here and 
    # use the center of the system and move things relative to that, but it turns out randomly 
    # moving in different directions works better. Oh well. That's what I get for trying to be cleve.r 
    get_fds=PDBParserObj.FD_coords
    # get loops
    fd_loops=PDBParserObj.FD_loop_coords
    for loop_num in fd_loops:
        get_fds[loop_num]=(fd_loops[loop_num][0])
    all_fds={}
    for fd in get_fds:
        fd_num=int(fd.split('_')[-1])
        all_fds[fd_num]=get_fds[fd]
    # sort
    fds_in_order=sorted(all_fds)

    # if linear placement is True...
    if linear_placement==True:
        PDBParserObj = center_fd_at_0(PDBParserObj)


    # only move fds that have IDRs between them. 
    num_fd_move=1
    for fd in fds_in_order:
        if fd+2 in fds_in_order:
            current_attempts=0
            success=False
            if verbose:
                print(f'Moving FDs {num_fd_move}')
                num_fd_move+=1
            # get coords for fds. 
            domain1_coords_ind = all_fds[fd]
            domain2_coords_ind = all_fds[fd+2]
            # list to restrict regions we check for clashing. Lets us only check clashing of
            # things we have already move, can ignore everything else. 
            check_coords=[]
            restrict_some_coords=False
            if len(all_fds)>2:
                for fd_ind in all_fds:
                    # only need to check things that are before what we are moving right now. 
                    if fd_ind <= fd:
                        check_coords.append(all_fds[fd_ind])
                        restrict_some_coords=True

            # make sure we have this empty nested list of we are checking everyhting
            if restrict_some_coords==False:
                check_coords=[[]]

            # Get actual coords for fd regions.
            domain1_coords=get_region_coords(PDBParserObj, specific_regions=[domain1_coords_ind], return_list=False)
            domain2_coords=get_region_coords(PDBParserObj, specific_regions=[domain2_coords_ind], return_list=False)            
            
            # for repopulating PDBParser after confirming everything is good..
            aa_and_at_dom2=[]
            domain2_coords_list=[]
            for aa in domain2_coords:
                for at in domain2_coords[aa]:
                    aa_and_at_dom2.append([aa, at])
                    domain2_coords_list.append(domain2_coords[aa][at])
            # get e2e objective.
            if run_prediction==False:
                objective_dist = len([aa for aa in range(domain1_coords_ind[-1], domain2_coords_ind[0])])*length_multiplier
            else:
                objective_dist = predict_e2e(PDBParserObj.sequence[domain1_coords_ind[-1]: domain2_coords_ind[0]])

            # make sure we update the CA coords as we go. 
            atoms=PDBParserObj.all_atom_coords_by_index
            # Convert input coordinates to NumPy arrays
            ref_point1 = np.array(atoms[domain1_coords_ind[-1]-1]['CA'])
            ref_point2 = np.array(atoms[domain2_coords_ind[0]]['CA'])  
            # while not at the max number...
            while current_attempts < total_attempts:
                current_attempts+=1
                if verbose==True:
                    print(f'on attempt {current_attempts}!')
                # get random point on sphere. This is the objective distance from the end of the first FD
                if linear_placement==False:
                    random_point = random_coordinate_on_sphere_surface(ref_point1, objective_dist)
                    # translate coords of fd2 by the distance we want in the random direction from FD1
                    translated_coords=translate_coordinates(domain2_coords_list, ref_point2, random_point)
                else:
                    # need to set translated coords at this part, probs gonna use a function
                    target_coord = (ref_point1[0]+objective_dist, ref_point1[1], ref_point1[2])
                    # translate the coords
                    translated_coords=translate_coordinates(domain2_coords_list, ref_point2, target_coord)
                    # only try this once because we are going to end up in the same location every time. 
                    current_attempts+=total_attempts

                # track numbers of aa and atom for repopulating dict
                raw_num_ind=0
                # temp copy of PDBParserObj to update and test for clashing
                temp_PDBParserObj=deepcopy(PDBParserObj)
                for coord in translated_coords:
                    # get aa num and atom name to repopulate thigns
                    cur_pair = aa_and_at_dom2[raw_num_ind]
                    cur_ats = temp_PDBParserObj.all_atom_coords_by_index[cur_pair[0]]
                    cur_ats[cur_pair[1]]=coord
                    # update all atom coords dict in pdbparserobj
                    temp_PDBParserObj.all_atom_coords_by_index[cur_pair[0]]=cur_ats
                    # update the index so we track atoms and aa numbers
                    raw_num_ind+=1

                # check for clashing
                if is_min_dist_satisfied_FDs(temp_PDBParserObj, domain2_coords_ind, 
                    check_these_coords=check_coords, CA_clash_dist=parameters.CA_clash_dist)!=True:
                        success=False
                else:
                    # if we succeeded, copy the PDBParserObj back from temp_PDBParserObj
                    PDBParserObj=deepcopy(temp_PDBParserObj)
                    success=True
                    break
    if success==False:
        raise dodoException('Unable to place the FD withut clashing!')
    # return updated obj with moved fds.
    return PDBParserObj 
            

def build_connecting_IDRs(PDBParserObj, bond_length=parameters.CA_bond_length,
    clash_dist=parameters.CA_clash_dist, attempts_per_coord=5000, verbose=False, debug=False):
    '''
    For buidling idrs after placing fds. Basically starts at an alpha carbon length
    from the first FD and builds to an alpha carbon length from the next FD. 

    Parameters
    ----------
    PDBParserObj : PDBParser
        PDBParser object
    bond_length : float
        length of bond between atoms. Default is those specified in parameters.py
    clash_dist : float
        minimum distance between CA atoms. Default is those specified in parameters.py
    attempts_per_coord : int
        number of attempts to try to place each coordinate. Default is 5000.
    verbose : bool
        if True, prints out why min distance was not satisfied. Default is False.
    debug : bool
        if True, prints out some extra stuff for debugging. Default is False.

    Returns
    -------
    PDBParserObj
        PDBParser object with IDRs built.
    '''
    # get all atoms by ind
    all_atoms = PDBParserObj.all_atom_coords_by_index
    # get fds
    get_fds=PDBParserObj.FD_coords
    # get loops
    fd_loops=PDBParserObj.FD_loop_coords
    for loop_num in fd_loops:
        get_fds[loop_num]=(fd_loops[loop_num][0])
    all_fds={}
    for fd in get_fds:
        fd_num=int(fd.split('_')[-1])
        all_fds[fd_num]=get_fds[fd]
    # sort
    fdnums=list(all_fds.keys())    
    idrs=PDBParserObj.IDR_coords
    idrnums=[]
    for idr in idrs:
        idrnums.append(int(idr.split('_')[-1]))

    # final radius distances
    final_distances={10: 15, 9: 13.8, 8: 12.5, 7: 11.2, 6: 10, 5: 8.8, 4: 7.5, 3: 6.2, 2: 5.1, 1: 3.8}
    off_by_dict={30: 10, 29: 9.5, 28: 9, 27: 8.5, 26: 8, 25: 8, 24: 8, 23: 8, 22: 8,21: 6.5,20:5, 19:5, 18:5,  17:5, 16:5, 
    15:5, 14:5, 13:5, 12: 5, 10: 4.5, 9: 4, 8: 4, 7: 3.5, 6: 2.5, 5: 2, 4: 1.5, 3: 1, 2: 1, 1: 0.75}
    # just idrs with +/- fds
    for idr_num in idrnums:
        if idr_num-1 in fdnums and idr_num+1 in fdnums:
            IDR_start_res_ind = all_fds[idr_num-1][1]
            IDR_end_res_ind = all_fds[idr_num+1][0]-1
            final_coordinate = all_atoms[IDR_end_res_ind+1]['CA']
            # get idr indices and IDR length
            all_idr_indices=[aa for aa in range(IDR_start_res_ind, IDR_end_res_ind+1)]
            idr_length=len(all_idr_indices)
            # get idr end to end dist
            IDR_end_to_end_dist = get_res_dist(all_atoms[IDR_start_res_ind-1]['CA'], all_atoms[IDR_end_res_ind+1]['CA'])
            # make radius distances
            if idr_length>10:
                max_dist_out=IDR_end_to_end_dist+((0.2)*(idr_length-10))
                radius_distances=sorted(np.linspace(14.8, max_dist_out, idr_length-10),reverse=True)
                # finish radius distances
                for i in range(10,0,-1):
                    radius_distances.append(final_distances[i])
            else:
                radius_distances=[]
                for i in range(idr_length,0,-1):
                    radius_distances.append(final_distances[i])
            radius_distances_by_res={}
            for ind_num, res_ind in enumerate(all_idr_indices):
                radius_distances_by_res[res_ind]=radius_distances[ind_num]            
            for idr_val in range(IDR_start_res_ind, IDR_end_res_ind+1):
                num_to_end=IDR_end_res_ind-idr_val+1
                success=False
                cur_attempt=0
                if num_to_end in list(off_by_dict.keys()):
                    off_by=off_by_dict[num_to_end]
                else:
                    off_by=None
                while cur_attempt < attempts_per_coord:
                    cur_attempt+=1
                    cur_radius_dist=radius_distances_by_res[idr_val]
                    prev_CA_coord=all_atoms[idr_val-1]['CA']
                    new_IDR_coord=random_walk_step(prev_CA_coord, bond_length)
                    # check for clashes
                    if is_min_distance_satisfied(PDBParserObj, idr_val, new_IDR_coord, 'CA', verbose=debug):
                        if is_point_within_sphere_constraints(final_coordinate, new_IDR_coord, cur_radius_dist, off_by=off_by, verbose=debug):
                            if debug:
                                print(get_res_dist(new_IDR_coord, final_coordinate))
                            success=True
                            #Update PDBParserObj with the new CA coord
                            PDBParserObj.all_atom_coords_by_index[idr_val]={'CA':new_IDR_coord}
                            # update all atoms
                            all_atoms=PDBParserObj.all_atom_coords_by_index
                            break
                if success==False:
                    raise dodoException(f'Unable to place the IDR without clashing! Failed on {idr_val}')
    return PDBParserObj

def build_terminal_IDRs(PDBParserObj, mode, verbose=True, 
    bond_length=parameters.CA_bond_length,clash_dist=parameters.CA_clash_dist, 
    attempts_per_coord=5000, debug=False):
    '''
    For building C and N terminal IDRs. Basically figures out either the predicted
    end to end dist of the IDR or multiplies the IDR length by a decimal value (depending
    on 'mode') and then places the starting or ending (depending on terminus) residue that 
    distance away from the FD such that the final end to end distance of the IDR is the 
    distance that we want it to be. 

    Parameters
    ----------
    PDBParserObj : PDBParser
        PDBParser object
    mode : string
        mode to use for placing FDs. Can be super_compact, compact, normal,
        expanded, super_expanded, max_expansion, or predicted.
    bond_length : float
        length of bond between atoms. Default is those specified in parameters.py
    clash_dist : float
        minimum distance between CA atoms. Default is those specified in parameters.py
    attempts_per_coord : int
        number of attempts to try to place each coordinate. Default is 5000.
    verbose : bool
        if True, prints out why min distance was not satisfied. Default is False.
    debug : bool
        if True, prints out some extra stuff for debugging. Default is False.

    Returns
    -------
    PDBParserObj
        PDBParser object with IDRs built.
    '''
    # check the mode
    run_prediction=False
    if mode =='super_compact':
        length_multiplier = 0.2
    elif mode == 'compact':
        length_multiplier=0.55
    elif mode == 'normal':
        length_multiplier=0.8
    elif mode == 'expanded':
        length_multiplier = 1.05
    elif mode == 'super_expanded':
        length_multiplier = 1.3
    elif mode == 'max_expansion':
        length_multiplier=1.65
    elif mode == 'predicted':
        run_prediction=True
    else:
        raise Exception('Invalid mode specified. Please specify super_compact, compact, normal, expanded, super_expanded, or max_expansion, or predicted.')

    # get all atoms by ind
    all_atoms = PDBParserObj.all_atom_coords_by_index

    # keep track of original attempts per coord
    original_attempts_per_coord=copy(attempts_per_coord)

    #get_all_regions
    all_regions=PDBParserObj.regions_dict
    terminal_regions=[list(all_regions.keys())[0],list(all_regions.keys())[-1]]
    for terminal_region in terminal_regions:
        if 'idr' in terminal_region:
            idr_coords=all_regions[terminal_region]
            idr_indices=[aa for aa in range(idr_coords[0], idr_coords[1])]
            if terminal_region==list(all_regions.keys())[0]:
                N_terminal_IDR=True
                adjust_res=1
                radius_center_residue_ind=idr_coords[-1]
            else:
                N_terminal_IDR=False
                adjust_res=-1
                radius_center_residue_ind=idr_coords[0]-1
                idr_indices.append(idr_indices[-1]+1)
            len_idr=len(idr_indices)
            # get e2e
            if run_prediction==False:
                IDR_end_to_end_dist=len_idr*length_multiplier
            else:
                IDR_end_to_end_dist=predict_e2e(PDBParserObj.sequence[idr_indices[0]: idr_indices[-1]+1])
            if N_terminal_IDR==False:
                idr_indices=sorted(idr_indices,reverse=True)
            if len_idr>10:
                    # final distances and off by vlaues for the sphere check
                final_distances={10: 15, 9: 13.8, 8: 12.5, 7: 11.2, 6: 10, 5: 8.8, 4: 7.5, 3: 6.2, 2: 5.1, 1: 3.8}
                off_by_dict={30: 10, 29: 9.5, 28: 9, 27: 8.5, 26: 8, 25: 8, 24: 8, 23: 8, 22: 8,21: 6.5,20:5, 19:5, 18:5,  17:5, 16:5, 
                            15:5, 14:5, 13:5, 12: 5, 10: 4.5, 9: 4, 8: 4, 7: 3.5, 6: 2.5, 5: 2, 4: 1.5, 3: 1, 2: 1, 1: 0.75}
                attempts_per_coord=original_attempts_per_coord
                max_dist_out=IDR_end_to_end_dist+((0.2)*(len_idr-10))
                radius_distances=sorted(np.linspace(14.8, max_dist_out, len_idr-10),reverse=True)
                # finish radius distances
                for i in range(10,0,-1):
                    radius_distances.append(final_distances[i])
            else:
                attempts_per_coord=20000
                radius_distances=[]
                # trying alternative attempt...
                dist_from_final = {10: 10.8, 9: 9.8, 8: 9.3, 7: 8.8, 6: 7.3, 5: 6.8, 4: 6.3, 3: 5.8, 2: 4.8, 1: 3.8}
                off_by_dict={10: 1.5, 9: 1.5, 8: 1.5, 7: 1.5, 6: 1.5, 5: 1.5, 4: 1.5, 3: 1, 2: 1, 1: 0.75}
                final_distances={10: 21.8, 9: 19.8, 8: 17.3, 7: 15.8, 6: 13.3, 5: 11.8, 4: 9.8, 3: 7.8, 2: 5.8, 1: 3.8}
                for i in range(len_idr,0,-1):
                    radius_distances.append(final_distances[i])
            radius_distances_by_res={}
            for ind_num, res_ind in enumerate(idr_indices):
                radius_distances_by_res[res_ind]=radius_distances[ind_num]  
            #try to place starting coord. 
            radius_center_coord=all_atoms[radius_center_residue_ind]['CA']
            for idr_index_num in idr_indices:
                success=False
                number_res = len_idr-idr_indices.index(idr_index_num)
                cur_radius_dist = radius_distances_by_res[idr_index_num]
                if number_res in list(off_by_dict.keys()):
                    off_by=off_by_dict[number_res]
                else:
                    off_by=None
                cur_attempt=0
                while cur_attempt < attempts_per_coord:
                    cur_attempt+=1
                    if idr_index_num==idr_indices[0]:
                        new_IDR_coord=random_coordinate_on_sphere_surface(radius_center_coord, IDR_end_to_end_dist)
                        check_sphere=False
                    else:
                        prev_CA_coord=all_atoms[idr_index_num-adjust_res]['CA']
                        new_IDR_coord=random_walk_step(prev_CA_coord, bond_length) 
                        check_sphere=True
                    if check_sphere==True:
                        within_sphere = is_point_within_sphere_constraints(radius_center_coord, new_IDR_coord, cur_radius_dist, off_by=off_by, verbose=debug)
                    else:
                        within_sphere=True
                    if len(idr_indices)<10:
                        if within_sphere==False:
                            new_IDR_coord = randomish_step(prev_CA_coord, radius_center_coord, dist_from_final[number_res],
                                                step_from_start_coord_dist=3.8)
                            within_sphere = is_point_within_sphere_constraints(radius_center_coord, new_IDR_coord, cur_radius_dist, off_by=off_by, verbose=debug)
                    if within_sphere==True:
                        if is_min_distance_satisfied(PDBParserObj, idr_index_num, new_IDR_coord, 'CA', verbose=debug):
                            PDBParserObj.all_atom_coords_by_index[idr_index_num]={'CA':new_IDR_coord}
                            # update all atoms
                            all_atoms=PDBParserObj.all_atom_coords_by_index
                            success=True
                            break
                    else:
                        success=False
                if success==False:
                    raise dodoException('Unable to place starting coordinate!')
                # have start coord. Time to go on!
    return PDBParserObj


def build_loops(PDBParserObj, verbose=True, 
        bond_length=parameters.CA_bond_length, 
        clash_dist=parameters.CA_clash_dist, 
        attempts_per_coord=10000, debug=False):
    '''
    Function for building loops in FDs. If there are multiple loops, it
    iterates through all the loops for that FD. 

    Parameters
    ----------
    PDBParserObj : PDBParser
        PDBParser object
    bond_length : float
        length of bond between atoms. Default is those specified in parameters.py
    clash_dist : float
        minimum distance between CA atoms. Default is those specified in parameters.py
    attempts_per_coord : int
        number of attempts to try to place each coordinate. Default is 5000.
    verbose : bool
        if True, prints out why min distance was not satisfied. Default is False.
    debug : bool
        if True, prints out some extra stuff for debugging. Default is False.

    Returns
    -------
    PDBParserObj
        PDBParser object with loops built.
    '''
    # final radius distances
    final_distances={10: 15, 9: 13.8, 8: 12.5, 7: 11.2, 6: 10, 5: 8.8, 4: 7.5, 3: 6.2, 2: 5.0, 1: 3.9}
    off_by_dict={53: 42.5, 52: 41, 51: 39.5, 50: 37, 49: 35.5, 48:34, 47: 32.5, 46:31, 45:29.5, 44:28,
    43: 26.5, 42:25, 41:23.5, 40: 22, 39: 20.5,38:20, 37:19.5, 36:19, 35:17.5,34: 16, 33:14.5, 32:13, 31:11.5, 30: 10, 
    29: 9.5, 28: 9, 27: 8.5, 26: 8, 25: 8, 24: 8, 23: 8, 22: 8,21: 6.5,20:5, 19:5, 18:5,  17:5, 16:5, 
    15:5, 14:5, 13:5, 12: 5, 10: 4.5, 9: 4, 8: 4, 7: 3.5, 6: 2.5, 5: 2, 4: 1.5, 3: 1, 2: 1, 1: 0.75}
    # get loops
    fd_loops=PDBParserObj.FD_loop_coords
    # iterate thorugh fds and loops. 
    for fd_loop_name in fd_loops:
        fd_and_loops=fd_loops[fd_loop_name]
        fd_region=fd_and_loops[0]
        cur_loops=fd_and_loops[1]
        # iterate through loops.
        for loop in cur_loops:
            loop_indices=[aa for aa in range(loop[0], loop[1]+1)]
            len_idr=len(loop_indices)
            if len_idr>10:
                min_dist = get_res_dist(PDBParserObj.all_atom_coords_by_index[loop_indices[0]-1]['CA'], PDBParserObj.all_atom_coords_by_index[loop_indices[-1]+1]['CA'])
                if (14.8+((len(loop_indices)-10)*0.8)) > min_dist:
                    radius_distances=sorted(np.linspace(14.8, (14.8+((len(loop_indices)-10)*0.8)), len_idr-10),reverse=True)
                    # finish radius distances
                    for i in range(10,0,-1):
                        radius_distances.append(final_distances[i])
                else:
                    radius_distances=sorted(np.linspace(3.8, min_dist+3, len_idr), reverse=True)
                
            else:
                radius_distances=[]
                for i in range(len_idr,0,-1):
                    radius_distances.append(final_distances[i])
            radius_distances_by_res={}
            for ind_num, res_ind in enumerate(loop_indices):
                radius_distances_by_res[res_ind]=radius_distances[ind_num]              
            # get loop indices.
            radius_center_coord=PDBParserObj.all_atom_coords_by_index[loop_indices[-1]+1]['CA']
            for resnum, loop_aa_ind in enumerate(loop_indices):
                success=False
                cur_radius_dist = radius_distances_by_res[loop_aa_ind]
                cur_radius_ind_num=len(loop_indices)-resnum
                cur_attempt=0
                if cur_radius_ind_num in list(off_by_dict.keys()):
                    off_by=off_by_dict[cur_radius_ind_num]
                else:
                    off_by=None
                while cur_attempt < attempts_per_coord:
                    cur_attempt+=1
                    prev_CA_coord=PDBParserObj.all_atom_coords_by_index[loop_aa_ind-1]['CA']
                    new_IDR_coord=random_walk_step(prev_CA_coord, bond_length)
                    # check for clashes
                    if is_min_distance_satisfied(PDBParserObj, loop_aa_ind, new_IDR_coord, 'CA', verbose=debug):
                        if is_point_within_sphere_constraints(radius_center_coord, new_IDR_coord, cur_radius_dist, off_by=off_by, verbose=debug):
                            success=True
                            #Update PDBParserObj with the new CA coord
                            PDBParserObj.all_atom_coords_by_index[loop_aa_ind]={'CA':new_IDR_coord}
                            # update all atoms
                            all_atoms=PDBParserObj.all_atom_coords_by_index
                            break
                if success==False:
                    raise dodoException(f'Unable to place the loop without clashing! Failed on {loop_aa_ind}')
    return PDBParserObj


def build_structure(PDBParserObj, mode='predicted', attempts_per_region=40, 
    attempts_per_coord=2500, linear_placement=False, verbose=True, very_verbose=False,
    num_models=1, beta_for_FD_IDR=False):
    '''
    Function for building the entire structure. Combines all the other functions
    in this module. 

    Parameters
    ----------
    PDBParserObj : PDBParser
        PDBParser object
    mode : string
        mode to use for placing FDs. Can be super_compact, compact, normal,
        expanded, super_expanded, max_expansion, or predicted.
    attempts_per_region : int
        number of attempts to try to place each FD. Default is 20.
    attempts_per_coord : int
        number of attempts to try to place each coordinate. Default is 2000.
    linear_placement : bool
        whether to place the folded domains across a linear axis. 
        Default : False
    verbose : bool
        if True, prints out why min distance was not satisfied. Default is False.
    very_verbose : bool
        if True, prints out a lot of stuff. Default is False. Mainly for debugging. 
    num_models : int
        number of models to build. This will rebuild the IDR multiple times
        but will hold the FDs constant. Will build the number of IDRs specified here.
        Default : 1
    beta_for_FD_IDR : bool
        whether to overwrite beta for the FD IDRs. FD = 0, IDR =100. Nifty for VMD.
        Default : False

    Returns
    -------
    PDBParserObj
        PDBParser object with FDs placed, IDRs connecting FDs built,
        any loops build, and terminal IDRs built.
    '''
    # dict to hold models:
    models = {}

    # remove the IDRs. This is so we don't need to worry about them 
    # clashing before we build everything else. 
    PDBParserObj=remove_IDRs_loops(PDBParserObj)
    # get all the regions we need to build
    all_regions=PDBParserObj.regions_dict
    # if we have FDs or FDs with loops...
    if len(list(PDBParserObj.FD_coords.keys()))+len(list(PDBParserObj.FD_loop_coords.keys()))>1:
        success=False
        cur_attempt=0
        if verbose==True:
            print('Setting position of folded domains')  
        # attempt to place the FDs.        
        while cur_attempt < attempts_per_region:
            cur_attempt+=1
            try:
                PDBParserObj=place_FDs(PDBParserObj,mode=mode, linear_placement=linear_placement,
                verbose=very_verbose, total_attempts=30)
                success=True
                break
            except dodoException:
                pass
        if success==False:
            raise dodoException('Unable to place FDs.')
    
    # make a deep copy of the starting PDBParserObj to use to overwrite stuff
    # if we have multiple models. 
    starting_PDBParserObj = deepcopy(PDBParserObj)

    for model in range(1, num_models+1):
        if model > 1:
            PDBParserObj = deepcopy(starting_PDBParserObj)

        # if we have loops, build in the loops.
        if PDBParserObj.FD_loop_coords!={}:
            success=False
            cur_attempt=0
            if verbose==True:
                print('Creating disordered loops.')        
            while cur_attempt < attempts_per_region:
                cur_attempt+=1
                try:
                    PDBParserObj=build_loops(PDBParserObj, verbose=very_verbose, attempts_per_coord=attempts_per_coord)
                    success=True
                    break
                except dodoException:
                    pass
            if success==False:
                raise dodoException('Unable to make loops.') 

        # if we have IDRs... *This doesn't need to worry if they are N or C terminal, 
        # because the function will figure that out. 
        if PDBParserObj.IDR_coords!={}:
            success=False
            cur_attempt=0
            if verbose==True:
                print('Connecting FDs with IDRs.')        
            while cur_attempt < attempts_per_region:
                cur_attempt+=1
                try:
                    PDBParserObj=build_connecting_IDRs(PDBParserObj, verbose=very_verbose, attempts_per_coord=attempts_per_coord)
                    success=True
                    break
                except dodoException:
                    pass
            if success==False:
                raise dodoException('Unable to generate IDRs to connect FDs.')            
        
        # if we have terminal IDRs, add those in.
        if 'idr' in list(all_regions.keys())[0] or 'idr' in list(all_regions.keys())[-1]:
            success=False
            cur_attempt=0
            if verbose==True:
                print('Adding N and / or C terminal IDRs.')        
            while cur_attempt < attempts_per_region:
                cur_attempt+=1
                try:
                    PDBParserObj=build_terminal_IDRs(PDBParserObj, mode=mode, verbose=very_verbose, attempts_per_coord=attempts_per_coord)
                    success=True
                    break
                except dodoException:
                    pass
            if success==False:
                raise dodoException('Unable to generate terminal IDRs.')          
        
        # get atom count. This is for CONECT lines later
        tot_atoms=0
        for aa_ind in PDBParserObj.all_atom_coords_by_index:
            tot_atoms+=len(PDBParserObj.all_atom_coords_by_index[aa_ind])

        # update atom count in PDBParserObj.
        PDBParserObj.number_atoms=tot_atoms
        # final check
        if len(PDBParserObj.all_atom_coords_by_index)!=len(PDBParserObj.sequence):
            raise dodoException('Length of coordinates by aa index does not match sequence length!')

        # get the sequence *using the all atom coordinates* so we can make sure nothing
        # got dropped as we were updating the coordinates using each function.
        built_sequence = get_seq_from_all_atom_coords(PDBParserObj)
        if built_sequence != PDBParserObj.sequence:
            raise dodoException('Built sequence does not match input sequence!')

        # if we are supposed to overwrite beta vals, do so
        if beta_for_FD_IDR==True:
            PDBParserObj=overwrite_beta_for_vis(PDBParserObj)
        
        # add to models dict
        models[model]=deepcopy(round_coordinates(PDBParserObj))
    

    # return coordinates rounded to 3 places.
    return models
    

def build_idr_from_sequence(sequence, mode, 
    end_coord=(0,0,0), attempts_per_idr=10,
    attempts_per_coord=5000,
    bond_length=parameters.CA_bond_length, 
    clash_dist=parameters.CA_clash_dist):
    '''
    Function for building an IDR from a sequence.

    sequence: string
        sequence of the IDR as a string
    mode : string
        mode to use for placing FDs. Can be super_compact, compact, normal,
        expanded, super_expanded, max_expansion, or predicted.
    end_coord : tuple
        the XYZ coordinate of the starting residue. Can be whatever you'd like..
    attempts_per_idr : int
        number of times to try to build the idr
    attempts_per_coord : int
        number of attempts to try to place each coordinate. Default is 5000.
    bond_length : float
        length of bond between atoms. Default is those specified in parameters.py
    clash_dist : float
        minimum distance between CA atoms. Default is those specified in parameters.py
    
    Returns
    -------
    dict
        dict with residues info and cooridnates. Basically info needed to
        make everything for the structure in build.py
    '''
    # check the mode
    run_prediction=False
    if mode =='super_compact':
        length_multiplier = 0.2
    elif mode == 'compact':
        length_multiplier=0.55
    elif mode == 'normal':
        length_multiplier=0.8
    elif mode == 'expanded':
        length_multiplier = 1.05
    elif mode == 'super_expanded':
        length_multiplier = 1.3
    elif mode == 'max_expansion':
        length_multiplier=1.65
    elif mode == 'predicted':
        run_prediction=True
    else:
        raise Exception('Invalid mode specified. Please specify super_compact, compact, normal, expanded, super_expanded, or max_expansion, or predicted.')
    
    # get the sequence length
    seq_len=len(sequence)

    # get the end to end distance
    if run_prediction==True:
        e2e=predict_e2e(sequence)
    else:
        e2e=seq_len*length_multiplier

    # final distances and off by vlaues for the sphere check
    final_distances={10: 15, 9: 13.8, 8: 12.5, 7: 11.2, 6: 10, 5: 8.8, 4: 7.5, 3: 6.2, 2: 5.0, 1: 3.9}
    off_by_dict={30: 10, 29: 9.5, 28: 9, 27: 8.5, 26: 8, 25: 8, 24: 8, 23: 8, 22: 8,21: 6.5,20:5, 19:5, 18:5,  17:5, 16:5, 
    15:5, 14:5, 13:5, 12: 5, 10: 4.5, 9: 4, 8: 4, 7: 3.5, 6: 2.5, 5: 2, 4: 1.5, 3: 1, 2: 1, 1: 0.75}

    # make radius distances
    if seq_len>10:
        max_dist_out=e2e+((0.2)*(seq_len-10))
        radius_distances=sorted(np.linspace(14.8, max_dist_out, seq_len-10),reverse=True)
        # finish radius distances
        for i in range(10,0,-1):
            radius_distances.append(final_distances[i])
    else:
        radius_distances=[]
        for i in range(seq_len,0,-1):
            radius_distances.append(final_distances[i])
    radius_distances_by_res={}
    for ind_num, res_ind in enumerate([aa for aa in range(0, len(sequence))]):
        radius_distances_by_res[res_ind]=radius_distances[ind_num]  
    
    # track if we successfully built the IDR.
    successful_build=False
    build_attempts=0

    while build_attempts < attempts_per_idr:
        # make list to hold coords and add first aa
        idr_coords=[random_coordinate_on_sphere_surface(end_coord, e2e)]

        # iterate through residues
        for res_ind in range(1,seq_len):
            cur_dist_ind=seq_len-res_ind
            prev_CA_coord=idr_coords[-1]
            success=False
            cur_attempt=0
            if cur_dist_ind in list(off_by_dict.keys()):
                off_by=off_by_dict[cur_dist_ind]
            else:
                off_by=None
            cur_radius_dist = radius_distances_by_res[res_ind]
            while cur_attempt < attempts_per_coord:
                new_IDR_coord=random_walk_step(prev_CA_coord, bond_length) 
                if is_point_within_sphere_constraints(end_coord, 
                    new_IDR_coord, cur_radius_dist, off_by=off_by, verbose=False):
                    if is_min_distance_satisfied_idr_from_seq(idr_coords, 
                        new_IDR_coord, clash_dist):
                        idr_coords.append(new_IDR_coord)
                        success=True
                        break
                cur_attempt+=1

        # if correct number of coords, we completed the build.
        if len(idr_coords)==seq_len:
            successful_build=True
            break
        # track number of build attempts
        build_attempts+=1

    # if we never succeeded, raise Exception.
    if successful_build==False:
        raise dodoException('Unable to build IDR!')
    
    # lists for stuffs
    res_names=[]
    res_indices=[]
    atom_indices=[]
    atom_names=[]
    CONECT_coords=[]

    if len(sequence)!=len(idr_coords):
        raise dodoException('The length of the sequence does not match the number of coordinates!')

    for ind, aa in enumerate(sequence):
        res_names.append(parameters.AADICT[aa])
        res_indices.append(ind+1)
        atom_indices.append(ind+1)
        atom_names.append('CA')
        if ind < len(sequence)-1:
            CONECT_coords.append([ind+1, ind+2])

    return {'xyz_list':idr_coords, 'residue_names':res_names, 'residue_indices':res_indices,
                'atom_indices':atom_indices, 'atom_names':atom_names, 'CONECT_coords':CONECT_coords}
