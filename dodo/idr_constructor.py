# DEPRECATED

# Bunch of functions here as I'm 
# trying different solutions for making IDRs that connect FDs.

import random
import math
import numpy as np
from sparrow import Protein as pr
from dodo.dodo_exceptions import dodoException

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


    theta = random.uniform(0, 2 * math.pi)
    phi = random.uniform(0, math.pi)

    # Calculate the Cartesian coordinates
    x = initial_coords[0] + distance * math.sin(phi) * math.cos(theta)
    y = initial_coords[1] + distance * math.sin(phi) * math.sin(theta)
    z = initial_coords[2] + distance * math.cos(phi)

    return (x, y, z)


def is_min_distance_satisfied(coords_list, coordinate_of_interest, min_distance=3.8):
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



def random_coordinate_with_biased_axis(start_coordinate, total_distance, 
    specified_distance, specified_coordinate):
    """
    Generate a random 3D coordinate that is a specific total distance from the starting coordinate
    and a specific distance in the X,Y, or Z direction from the starting coordinate.

    This isn't currently being implemented, but I'm keeping it here in case down the line I want
    to implement it for trickier structures.

    Paramters
    ---------
    start_coordinate : tuple
        The starting 3D coordinate (x0, y0, z0).
    total_distance : float
        The total distance between the starting coordinate and the random coordinate.
    specified_distance : float
        The distance in the X, Y, or Z direction between the starting coordinate and the random coordinate.
    specified_coordinate : str
        The coordinate (X, Y, or Z) in which the specified distance is measured.

    Returns:
        tuple: A pseudo-random 3D coordinate (x, y, z). Pseudo-random in that we can 
        bias it in any single direction by any amount less than 'total_distance'.
    """
    if len(start_coordinate) != 3:
        raise ValueError("Starting coordinate must be 3D, i.e., have three components (x0, y0, z0).")

    if specified_distance > total_distance:
        raise ValueError("X distance cannot be greater than the total distance.")

    # Calculate the remaining distance in the Y and Z directions
    remaining_distance = math.sqrt(total_distance**2 - specified_distance**2)

    # Generate random values for Y and Z coordinates within the sphere
    second_unit = random.uniform(-1, 1)
    third_unit = random.uniform(-1, 1)

    # Normalize the random values and scale them by the remaining distance
    normalization_factor = remaining_distance / math.sqrt(second_unit**2 + third_unit**2)
    second_unit *= normalization_factor
    third_unit *= normalization_factor

    # Calculate the final coordinates. What I did here was basically use a list
    # to hold required coordinates, remove the specified on, and then randomly assign
    # the remaining two coordinates to second_unit or third_unit. Keeps things 
    # agnostic as far as what coordinate we specify to bias. Yes, there's less code
    # that could accomplish this, but I made this whole thing in like 2 days so... sorry about that!
    required_coordinates = ['X', 'Y', 'Z']
    required_coordinates.remove(specified_coordinate)
    plane_to_loc = {'X':0, 'Y':1, 'Z':2}
    second_coord=random.choice(required_coordinates)
    required_coordinates.remove(second_coord)
    third_coord=required_coordinates[0]

    final_coords=[]
    for i in range(3):
        if i==plane_to_loc[specified_coordinate]:
            final_coords.append(start_coordinate[i]+specified_distance)
        elif i==plane_to_loc[second_coord]:
            final_coords.append(start_coordinate[i]+second_unit)
        elif i==plane_to_loc[third_coord]:
            final_coords.append(start_coordinate[i]+third_unit)
        else:
            raise dodoException('Problem with random specified coordinate fucntion.')

    return final_coords

def is_point_inside_sphere(point, center, radius):
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
    return distance <= radius

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

def generate_coordinates_on_line(coord1, coord2, distance):
    """
    Generate a 3D coordinates along a straight line, given two coordinates 
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
    direction_vector = coord2 - coord1

    # Calculate the unit vector in the direction of the line
    unit_vector = direction_vector / np.linalg.norm(direction_vector)

    # Generate coordinate along the line at specified interval
    return coord1+(distance*unit_vector)



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

    Returns:
        tuple: The center coordinate (x_center, y_center, z_center).
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


def find_closest_coordinate(coordinates, specific_coordinate):
    """
    Find the closest coordinate from a list of 3D coordinates to a specific coordinate.

    Parameters
    ----------:
    coordinates : list of tuples
        List of 3D coordinates [(x1, y1, z1), (x2, y2, z2), ...].
    specific_coordinate : tuple 
        The specific 3D coordinate to which distances are calculated.

    Returns
    --------
    tuple : 
        The closest coordinate from the list.
    """
    if not coordinates:
        raise ValueError("List of coordinates is empty.")

    if len(specific_coordinate) != 3:
        raise ValueError("The specific coordinate must be a 3D coordinate (x, y, z).")

    closest_coordinate = None
    closest_distance = float('inf')

    for coord in coordinates:
        if len(coord) != 3:
            raise ValueError("All coordinates in the list must be 3D, i.e., have three components (x, y, z).")

        distance = math.sqrt((coord[0] - specific_coordinate[0])**2 + (coord[1] - specific_coordinate[1])**2 + (coord[2] - specific_coordinate[2])**2)

        if distance < closest_distance:
            closest_distance = distance
            closest_coordinate = coord

    return closest_coordinate

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


def build_idr_between_fds(all_structure_coords, fd1_coords, fd2_coords, idr_seq, 
    mode='predicted', attempts_per_residue=10000, bond_length=3.8,
    clash_dist=3, custom_dist=None):
    '''
    Function to build an IDR that connects two FDs. 

    Parameters
    ----------
    all_structure_coords : list
        all of the coordinates.
    fd1_coords : list
        list of coordinates for first folded domain
    fd2_coords : list
        list of coordinates for second folded domain
    idr_seq : string
        the amino acid sequence of the IDR as a string
    mode : string
        the mode to use for generating the IDR.
        options are predicted, extended, or custom.
        predicted : uses the predicted end to end distance
        extended : uses the length of the IDR * 1.2
        compact : uses the length of the IDR * 0.5
        custom : uses the length of the IDR * custom_dist
    attempts_per_residue : int
        the number of attempts to generate a coordinate for each residue
    bond_length : float
        the bond length to use for the random walk
    clash_dist : float
        the minimum distance to keep between atoms
    custom_dist : float
        the custom distance to use if mode is custom.

    Returns
    -------
    final_coords : list
        A list that is the coords (in order) FD1, 
        the IDR coords generated here, and then FD2.

    '''
    # distance between FDs which will be how expanded the IDR is.
    if mode=='predicted':
        dist=predict_e2e(idr_seq)
    elif mode == 'extended':
        dist=len(idr_seq)*(1.2)
    elif mode =='compact':
        dist=len(idr_seq)*0.5
    elif mode=='custom':
        if custom_dist==None:
            raise Exception('Please specify custom_dist parameter if using custom mode!')
        dist=len(idr_seq)*custom_dist
    else:
        raise Exception('Options for mode are predicted, extended, compact, or custom!')

    # figure out if it's even theoretically possible to build the IDR.
    if dist > (len(idr_seq)*3.8): 
        raise dodoException('Cannot build IDR. Distance between FDs is too large.')

    # listof all coords. Use this to make sure that we don't 
    # build on top of previously generated coords.
    all_coords=[]
    all_coords.extend(fd1_coords)
    all_coords.extend(fd2_coords)

    # list for final coords. Using this to build the final returned coordinates
    final_coords=[]
    final_coords.extend(fd1_coords)

    # get the coordinates for approximate center of each fd.
    # this is to make sure we move the FDs correctly relative to each other.
    center_fd1=calculate_center_coordinate(fd1_coords)
    center_fd2=calculate_center_coordinate(fd2_coords)

    #get closest coords to centers.
    coord_fd1 = find_closest_coordinate(fd1_coords, center_fd1)
    coord_fd2=find_closest_coordinate(fd2_coords, center_fd2)

    # Figure out coordinate to place the last residue of FD2
    # on such that it is the appropriate distance from FD1.
    start_coord_fd2=generate_coordinates_on_line(coord_fd1, coord_fd2, dist)

    # add rest of FD2 for moving everything.....
    post_fd2_coords = all_structure_coords[(all_structure_coords.index(fd2_coords[-1])+1):]
    fd2_coords.extend(post_fd2_coords)

    # translate fd2 to new location along unit vector 
    # connectong approximate center of FD1 and FD2
    moved_fd2=translate_coordinates(fd2_coords, fd2_coords[0], start_coord_fd2)

    # Dict of final distances. Basically this lets us specify the sphere radius
    # for the final steps to build the IDR. I determined this more or less empirically
    # and found that without doing this, building the IDR is only successful ~5% of the time.
    final_distances={1:3.8, 2:5.4, 3:6.8, 4:8, 5:8.8, 6:9.6, 7:10.2, 8:10.6, 9:10.8, 10:11}
    if len(idr_seq)>10:
        radius_distances=sorted(np.linspace(12, bond_length+dist, len(idr_seq)-10),reverse=True)
        # finish radius distances
        for i in range(10,0,-1):
            radius_distances.append(final_distances[i])
    else:
        radius_distances=[]
        for i in range(len(idr_seq),0,-1):
            radius_distances.append(final_distances[i])

    
    # Now need to generate locations for each amino acid of the IDR.
    # list to hold IDR coords. Let's us hold something we can use as
    # our starting coord as we build the IDR.
    IDR_coords=[]

    # for each res,
    for res in range(0,len(idr_seq)):
        # reset attempt for current res
        cur_attempt=0

        # if first res of IDR, start coordinate is last res
        # of first FD. Otherwise, is last res of IDR.
        if res==0:
            start_coordinate=fd1_coords[-1]
        else:
            start_coordinate=IDR_coords[-1]

        # keep track of if we managed to do it!
        success=False
        # figure out radius size for this res. This is the radius 
        # of the sphere starting at moved_fd2[0], ie. the first residue
        # of the second FD after we translated it relative to FD1.
        radius_size=radius_distances[res]
        # while our number of attempts is under the max attempts per residue
        while cur_attempt < attempts_per_residue:
            # add 1 to current attempts
            cur_attempt+=1
            # take a random step.
            cur_step = random_walk_step(start_coordinate, bond_length)
            # check if within sphere.
            if is_point_inside_sphere(cur_step, moved_fd2[0], radius_size)==True:
                # check for any clashes
                if is_min_distance_satisfied(all_coords, cur_step, min_distance=clash_dist)==True:
                    # add to IDR coords and all coords. We don't really need 2 lists but I wrote
                    # this pretty quickly so I'm rolling with it.
                    IDR_coords.append(cur_step)
                    all_coords.append(cur_step)
                    # mark this as successful. 
                    success=True
                    break
        # If we get to any residue in the IDR and do not
        # succeed, we will raise an exception. 
        if success==False:
            raise dodoException('Could not generate a step without collisions after attempts.')
    # now add IDR_coords and then moved_fd2 coords to final coordinates.
    final_coords.extend(IDR_coords)
    final_coords.extend(moved_fd2)
    return final_coords
    

def build_idr(fd_coords, IDR_relative_terminus, idr_seq, 
    mode='predicted', attempts_per_residue=40000, bond_length=3.8,
    clash_dist=3, custom_dist=None):
    '''
    Function to build an IDR that does not connect FDs (N or C terminal IDRs)

    Parameters
    ----------
    fd_coords : list
        list of coordinates for folded domain
    IDR_relative_terminus : string
        the relative position of the IDR to the FD. If the IDR is amino acids 
        0-100 and the FD is 101-200, use 'N'. If the opposite (the IDR is after
        the FD), use 'C'
    idr_seq : string
        the amino acid sequence of the IDR as a string
    mode : string
        the mode to use for generating the IDR.
        options are predicted, extended, or custom.
        predicted : uses the predicted end to end distance
        extended : uses the length of the IDR * 1.2
        compact : uses the length of the IDR * 0.5
        custom : uses the length of the IDR * custom_dist
    attempts_per_residue : int
        the number of attempts to generate a coordinate for each residue
    bond_length : float
        the bond length to use for the random walk
    clash_dist : float
        the minimum distance to keep between atoms
    custom_dist : float
        the custom distance to use if mode is custom.

    Returns
    -------
    final_coords : list
        A list that is the coords (in order) FD1, 
        the IDR coords generated here, and then FD2.

    '''
    # distance between FDs which will be how expanded the IDR is.
    if mode=='predicted':
        dist=predict_e2e(idr_seq)
    elif mode == 'extended':
        dist=len(idr_seq)*(1.2)
    elif mode =='compact':
        dist=len(idr_seq)*0.5
    elif mode=='custom':
        if custom_dist==None:
            raise Exception('Please specify custom_dist parameter if using custom mode!')
        dist=len(idr_seq)*custom_dist
    else:
        raise Exception('Options for mode are predicted, extended, compact, or custom!')

    if IDR_relative_terminus=='N':
        start_coordinate = fd_coords[0]
        # original start coord
        original_start = fd_coords[0]
    elif IDR_relative_terminus=='C':
        start_coordinate = fd_coords[-1]
        # original start coord because 'start_coordinate' changes as we build the IDR.
        original_start = fd_coords[-1]
    else:
        raise Exception('Please specify whether IDR_relative_terminus is N or C.')

    # make radius distances...
    radius_distances=[3.8, 5.4, 6.8, 8, 8.8, 9.6, 10.2, 10.6, 10.8, 11]
    if len(idr_seq)>10:
        radius_distances.extend(np.linspace(12, bond_length+dist, len(idr_seq)-10))
    else:
        radius_distances=radius_distances[:len(idr_seq)]
        
    # track all coords for clashes.
    all_coords=[]
    all_coords.extend(fd_coords)

    # keep track of IDR coords specifically so we can reference as we build
    IDR_coords=[]
    # for residues in the IDR
    for res in range(0,len(idr_seq)):
        # reset attempt count
        cur_attempt=0
        # if not the first res of the IDR, the start coordinate is the
        # last residue we added to the IDR.
        if res!=0:
            start_coordinate=IDR_coords[-1]
        # get radius size. This lets us choose how expanded the IDR is.
        radius_size = radius_distances[res]
        # keep track of if we managed to do it!
        success=False
        # while we haven't attempted more than the max number of times..
        while cur_attempt < attempts_per_residue:
            # track attempt umber
            cur_attempt+=1
            # take a random step.
            cur_step = random_walk_step(start_coordinate, bond_length)
            # check if within sphere. As we build the IDR, the sphere gets bigger relative
            # to the start point. Let's us control how compact it is. 
            if is_point_inside_sphere(cur_step, original_start, radius_size)==True:
                # check for clashes
                if is_min_distance_satisfied(all_coords, cur_step, min_distance=clash_dist)==True:
                    # append to IDR_coords and all_coords
                    IDR_coords.append(cur_step)
                    all_coords.append(cur_step)
                    success=True
                    break
        # if we failed, raise exception.
        if success==False:
            raise dodoException('Could not generate a step without collisions after attempts.')
    # add the bits together. 
    final_coords=[]
    if IDR_relative_terminus=='N':
        # reverse idr coords
        reversed_coords=[]
        for coord in range(len(IDR_coords)-1,-1,-1):
            reversed_coords.append(IDR_coords[coord])
        final_coords.extend(reversed_coords)
        final_coords.extend(fd_coords)
    else:
        final_coords.extend(fd_coords)
        final_coords.extend(IDR_coords)
    return final_coords


def build_loop(fd_and_loop_xyz, loop_coords,
    attempts_per_residue=50000, bond_length=3.8, clash_dist=3):
    '''
    Function to build an IDR that connects two FDs. 

    Parameters
    ----------
    fd_and_loop_xyz : list
        list of coordinates for the entire folded domain including the loops.
    loop_coords : list
        The coordinates specifically for the loop within the FD.
    attempts_per_residue : int
        the number of attempts to generate a coordinate for each residue
    bond_length : float
        the bond length to use for the random walk
    clash_dist : float
        the minimum distance to keep between atoms

    Returns
    -------
    final_coords : list
        A list that is the coords (in order) FD1, 
        the IDR coords generated here, and then FD2.

    '''
    # figure out if it's even theoretically possible to build the loop.
    last_res_before_loop = fd_and_loop_xyz[loop_coords[0]-1]
    first_res_after_loop = fd_and_loop_xyz[loop_coords[-1]]
    min_theoretical_dist=get_res_dist(last_res_before_loop, first_res_after_loop)
    num_res_in_loop=len([aa for aa in range(loop_coords[0], loop_coords[1])])

    if min_theoretical_dist > (num_res_in_loop*3.8):
        raise dodoException('Cannot build loop. Distance between start of loop and end of loop is too large.')

    # break apart FD.
    fd_part_1 = fd_and_loop_xyz[:loop_coords[0]]
    fd_part_2 = fd_and_loop_xyz[loop_coords[1]:]
    
    # hold all cords except for loop.
    all_coords=[]
    all_coords.extend(fd_part_1)
    all_coords.extend(fd_part_2)

    # list for final coords. Using this to build the final returned coordinates
    final_coords=[]
    final_coords.extend(fd_part_1)

    # Dict of final distances. Basically this lets us specify the sphere radius
    # for the final steps to build the IDR. I determined this more or less empirically
    # and found that without doing this, building the IDR is only successful ~5% of the time.
    final_distances={1:3.8, 2:5.8, 3:7.8, 4:9.8, 5:11.8, 6:13.8, 7:15.8, 8:17.8, 9:19.8, 10:21.8}
    if num_res_in_loop>10:
        radius_distances=sorted(np.linspace(22, (22+(bond_length*(num_res_in_loop-10))), num_res_in_loop-10),reverse=True)
        # finish radius distances
        for i in range(10,0,-1):
            radius_distances.append(final_distances[i])
    else:
        radius_distances=[]
        for i in range(len(idr_seq),0,-1):
            radius_distances.append(final_distances[i])
    
    # Now need to generate locations for each amino acid of the loop.
    # list to hold loop coords. Let's us hold something we can use as
    # our starting coord as we build the IDR.
    loop_coords=[]

    # for each res,
    for res in range(0,num_res_in_loop):
        # reset attempt for current res
        cur_attempt=0

        # if first res of IDR, start coordinate is last res
        # of first FD. Otherwise, is last res of IDR.
        if res==0:
            start_coordinate=fd_part_1[-1]
        else:
            start_coordinate=loop_coords[-1]

        # keep track of if we managed to do it!
        success=False
        # figure out radius size for this res. This is the radius 
        # of the sphere starting at moved_fd2[0], ie. the first residue
        # of the second FD after we translated it relative to FD1.
        radius_size=radius_distances[res]


        # while our number of attempts is under the max attempts per residue
        while cur_attempt < attempts_per_residue:
            # add 1 to current attempts
            cur_attempt+=1
            # take a random step.
            cur_step = random_walk_step(start_coordinate, bond_length)
            # check if within sphere.
            if is_point_inside_sphere(cur_step, fd_part_2[0], radius_size)==True:
                # check for any clashes
                if is_min_distance_satisfied(all_coords, cur_step, min_distance=clash_dist)==True:
                    # add to IDR coords and all coords. We don't really need 2 lists but I wrote
                    # this pretty quickly so I'm rolling with it.
                    loop_coords.append(cur_step)
                    all_coords.append(cur_step)
                    # mark this as successful. 
                    success=True
                    break
        # If we get to any residue in the IDR and do not
        # succeed, we will raise an exception. 
        if success==False:
            raise dodoException('Could not generate a step without collisions after attempts.')
    # now add IDR_coords and then moved_fd2 coords to final coordinates.
    final_coords.extend(loop_coords)
    final_coords.extend(fd_part_2)
    return final_coords
    
def bond_length_check(coords, cutoff=4.2):
    '''
    function to check bond length between each amino acid
    across a set of coordinates.
    '''
    for i in range(0,len(coords)-1):
        if get_res_dist(coords[i],coords[i+1])>cutoff:
            return False
    return True


def build_new_structure(pdb_file_dict, regions_dict, IDR_connecting_FD_mode = 'predicted',
    IDR_connecting_FD_mode_custom=None, IDR_mode='predicted', IDR_custom = None,
    per_residue_attempts={'IDR':100000, 'loop':100000, 'connecting_IDR':100000}, 
    bond_length=3.8, clash_dist=3):
    '''
    What we've all been waiting for - the final function to build the structure.
    Boy howdy did I write this inefficently. But it works so good enough for now.
    Main reason to separate this all out as it goes is to avoid having to check all 
    coordinates for clashes every time we add something but... yikes.

    Parameters
    ----------
    pdb_file_dict : dict
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

    regions_dict : dict
        a dictionary of the regions to build.
            keys are a bit dynamic in that they are dependent on the structure and
            what regions it has. Basically, anything with an IDR should have
            the key be IDR_# and value is a list with the coords defining the IDR
            ex. {'idr_1':[0:20]}
            for a loop, key should be fd_with_loop_# and the value should be a nested
            list where the first element is the full FD + loop and the second is a list
            of lists containing coordinates for any loops. 
            ex. {'fd_with_loop_2' : [[100:200], [[110:120], [150,160]]]}
            if it's an FD, should be named with a key that has 'folded' in the name 
            and a value that is a list defining its location. 
            ex. {'folded_5':[960:1000]}

    IDR_connecting_FD_mode : string
        the mode to use for generating the connecting IDR.
        options are predicted, extended, compact, or custom.
        predicted : uses the predicted end to end distance
        extended : uses the length of the IDR * 1.2
        compact : uses the length of the IDR * 0.5
        custom : uses the length of the IDR * custom_dist

    IDR_connecting_FD_mode_custom : float
        the custom distance to use if IDR_connecting_FD_mode is custom.

    IDR_mode : string
        the mode to use for generating the IDRs.
        options are predicted, extended, compact, or custom.
        predicted : uses the predicted end to end distance
        extended : uses the length of the IDR * 1.2
        compact : uses the length of the IDR * 0.5
        custom : uses the length of the IDR * custom_dist

    IDR_custom : float
        the custom distance to use if IDR_mode is custom.

    per_residue_attempts : dict
        a dictionary of the number of attempts to generate a coordinate for each residue
        for each region type. 
        ex. {'IDR':40000, 'loop':50000, 'connecting_IDR':10000}

    bond_length : float
        the bond length to use for the random walk

    clash_dist : float
        the minimum distance to keep between atoms
    '''
    # first need to figure out which regions to build. Prioritize first building loops.
    all_regions=regions_dict.keys()
    loops=[region for region in all_regions if 'loop' in region]
    # get seq len
    seq_len = len(pdb_file_dict['sequence'])
    # get the CA coords
    ca_coords=pdb_file_dict['ca_coords']

    # build in the loops. This only checks clashes with the loop and the FD
    # that has the loop
    for loop in loops:
        current_set_of_loops = regions_dict[loop][1]
        current_fd = regions_dict[loop][0]
        current_fd_coords = ca_coords[current_fd[0]:current_fd[1]]
        for loop in current_set_of_loops:
            loop_added = build_loop(current_fd_coords, loop, 
                attempts_per_residue=per_residue_attempts['loop'],
                bond_length=bond_length, clash_dist=clash_dist)
            # over write coords for each sub loop in current set.
            current_fd_coords=loop_added
        # over write ca coords with loops added.
        ca_coords = ca_coords[:current_fd[0]] + current_fd_coords + ca_coords[current_fd[1]:]
        if len(ca_coords) != seq_len:
            raise dodoException('Wrong sequence length! Something went wrong with adding loops!')
        if bond_length_check(ca_coords)==False:
            raise dodoException('Wrong bond length! Something went wrong with adding loops!')

    # now build the IDRs connecting any FDs separated by IDRs. 
    # to keep things fast, this only checks the clashes between the 
    # two FDs and the IDR being bridged between them
    # get all IDRs
    IDRs=[region for region in all_regions if 'idr' in region]
    terminal_IDRs={}
    # for each IDR
    for IDR in IDRs:
        # if terminal N-terminal or C-terminal IDR
        if IDR == list(all_regions)[0]:
            terminal_IDRs['N']=IDR
        elif IDR == list(all_regions)[-1]:
            terminal_IDRs['C']=IDR
        else:
            # if now, gotta connect the FD.
            IDR_loc = list(regions_dict).index(IDR)
            
            preceeding_FD_regions = regions_dict[list(regions_dict)[IDR_loc-1]]
            following_FD_regions = regions_dict[list(regions_dict)[IDR_loc+1]]

            if type(preceeding_FD_regions[0])==list:
                preceeding_FD_regions=preceeding_FD_regions[0]
            if type(following_FD_regions[0])==list:
                following_FD_regions=following_FD_regions[0]

            cur_IDR_reg=regions_dict[IDR]

            current_IDR_seq = pdb_file_dict['sequence'][cur_IDR_reg[0]:cur_IDR_reg[1]]
            fd_1_coords = ca_coords[preceeding_FD_regions[0]:preceeding_FD_regions[1]]
            fd_2_coords = ca_coords[following_FD_regions[0]:following_FD_regions[1]]

            IDR_added = build_idr_between_fds(ca_coords, fd_1_coords, fd_2_coords, 
                current_IDR_seq, mode=IDR_connecting_FD_mode, 
                attempts_per_residue=per_residue_attempts['connecting_IDR'],
                bond_length=bond_length, clash_dist=clash_dist,
                custom_dist=IDR_connecting_FD_mode_custom)
            # over write ca coords with IDR added.
            ca_coords=ca_coords[:preceeding_FD_regions[0]] + IDR_added
            if len(ca_coords) != seq_len:
                raise dodoException('Wrong sequence length! Something went wrong with adding IDRs between FDs!')
            if bond_length_check(ca_coords)==False:
                raise dodoException('Wrong bond length! Something went wrong with adding IDRs between FDs!')


    # now build terminal IDRs. Because there's only 2 at most (by definitation),
    # we now will check for all clashes as we build these. By doing things in this
    # objectively obnoxious order, we keep this all working fast at the expense of
    # significant code bloat... but oh well.
    if 'N' in list(terminal_IDRs.keys()):
        IDR_coords = regions_dict[terminal_IDRs['N']]
        IDR_seq = pdb_file_dict['sequence'][IDR_coords[0]:IDR_coords[1]]
        fd_coords = ca_coords[IDR_coords[1]:]
        
        IDR_added = build_idr(fd_coords, 'N', IDR_seq, mode=IDR_mode, 
            attempts_per_residue=per_residue_attempts['IDR'],
            bond_length=bond_length, clash_dist=clash_dist,
            custom_dist=IDR_custom)
        # over write ca coords with IDR added.
        ca_coords = IDR_added
        if len(ca_coords) != seq_len:
            raise dodoException('Wrong sequence length! Something went wrong with adding N-terminal IDR!')
        if bond_length_check(ca_coords)==False:
            raise dodoException('Wrong bond length! Something went wrong with adding N-terminal IDR!')
        
    if 'C' in list(terminal_IDRs.keys()):
        IDR_coords = regions_dict[terminal_IDRs['C']]
        IDR_seq = pdb_file_dict['sequence'][IDR_coords[0]:IDR_coords[1]+1]
        fd_coords = ca_coords[:IDR_coords[0]]
        IDR_added = build_idr(fd_coords, 'C', IDR_seq, mode=IDR_mode, 
            attempts_per_residue=per_residue_attempts['IDR'],
            bond_length=bond_length, clash_dist=clash_dist,
            custom_dist=IDR_custom)
        # over write ca coords with IDR added.
        ca_coords = IDR_added
        if len(ca_coords) != seq_len:
            raise dodoException('Wrong sequence length! Something went wrong with adding C-terminal IDR!')
        if bond_length_check(ca_coords)==False:
            raise dodoException('Wrong bond length! Something went wrong with adding C-terminal IDr!')

    # return the final product
    if len(ca_coords) != seq_len:
        raise dodoException('Wrong sequence length! Something went wrong with adding IDRs between FDs!')
    if bond_length_check(ca_coords)==False:
        raise dodoException('Wrong bond length! Something went wrong with adding IDRs between FDs!')
    return ca_coords

        