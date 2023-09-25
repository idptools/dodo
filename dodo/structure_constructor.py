# rewrite of IDR_constructor.py to be a bit... better. 
# Noticed some issues in stability so going to try to fix that on a rewrite. 


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

def generate_coordinate_on_line(coord1, coord2, distance):
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

def random_coordinate_on_sphere_surface(center_coordinate, radius):
    """
    Generate a random 3D coordinate on the surface of a sphere with a specified radius.

    Args:
        center_coordinate (tuple): The center of the sphere (x_center, y_center, z_center).
        radius (float): The radius of the sphere.

    Returns:
        tuple: A random 3D coordinate on the sphere's surface.
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

def bond_length_check(coords, min_dist=2.8, max_dist=4.4):
    '''
    function to check bond length between each amino acid
    across a set of coordinates. Length must be within min dist and max dist.
    '''
    for i in range(0,len(coords)-1):
        if get_res_dist(coords[i],coords[i+1])>max_dist or get_res_dist(coords[i],coords[i+1])<min_dist:
            return False
    return True

def build_idr_from_sequence(IDR_end_to_end_dist, sequence,
    bond_length=3.8, clash_dist=3.4, attempts_start_coord=10000, 
    attempts_all_coords=30000, min_bond_dist = 2.8, max_bond_dist=4.4, start_coord=[0,0,0]):
    # just an IDR sequence and then the resultant IDR.
    if len(start_coord)!=3:
        raise dodoException('Start coordinate must be a 3D coordinate.')

    # end coordinate is a random distance away from the start coord based on the radius.     
    end_coord = random_coordinate_on_sphere_surface(start_coord, IDR_end_to_end_dist)

    IDR_length=len(sequence)

    # a few quick checks. 
    if IDR_end_to_end_dist > IDR_length*bond_length:
        raise dodoException('IDR length is not sufficient to be specified IDR_end_to_end_dist!')

    # list to hold IDR coorindates
    IDR_coords=[]
    IDR_coords.append(start_coord)

    # Dict of final distances. Basically this lets us specify the sphere radius
    # for the final steps to build the IDR. I determined this more or less empirically
    # and found that without doing this, building the IDR is only successful ~5% of the time.
    final_distances={10: 13, 9: 11.8, 8: 11.6, 7: 11.2, 6: 10.6, 5: 9.8, 4: 9, 3: 8, 2: 6, 1: 3.8}
    if IDR_length > 10:
        radius_distances=sorted(np.linspace(14, IDR_end_to_end_dist, IDR_length-10),reverse=True)
        # finish radius distances
        for i in range(10,0,-1):
            radius_distances.append(final_distances[i])
    else:
        radius_distances=[]
        for i in range(IDR_length,0,-1):
            radius_distances.append(final_distances[i])    

    # dict for distance per residue.
    dist_per_residue = {}
    for i in range(0, IDR_length):
        dist_per_residue[i]=radius_distances[i]
    
    # now can start from start coord and go from there with ever decreasing sphere size.
    for cur_IDR_res in range(1, IDR_length):
        current_radius_dist = dist_per_residue[cur_IDR_res]
        current_start = IDR_coords[-1]
        current_attempt=0
        # keep track of we succeeded. 
        success=False
        for attempt in range(0, attempts_all_coords):
            current_attempt+=1
            new_IDR_coord=random_walk_step(current_start, bond_length)
            if is_min_distance_satisfied(IDR_coords, new_IDR_coord, clash_dist):
                if is_point_inside_sphere(new_IDR_coord, end_coord, current_radius_dist):
                    IDR_coords.append(new_IDR_coord)
                    success=True
                    break
        if success==False:
            raise dodoException(f'Unable to build IDR residue {cur_IDR_res} after {current_attempt} attempts.')

    # final check for bond lengths.
    if bond_length_check(IDR_coords, min_dist=min_bond_dist, 
        max_dist=max_bond_dist)==False:
        raise dodoException('Error! bond_length_check function detected a bond length outside of the min and max bond length range!')

    # verify length is good
    if len(IDR_coords) != len(sequence):
        raise dodoException('Error! Length of final coordinates not equal to length of input sequence')

    # overwrite final coordinates in pdb dict
    return IDR_coords
    

def build_IDR_specific_dist_no_FDs(IDR_end_to_end_dist, pdb_dict,
    bond_length=3.8, clash_dist=3.4, attempts_start_coord=10000, 
    attempts_all_coords=30000, min_bond_dist = 2.8, max_bond_dist=4.4):
    # for building IDRs where there are no FDs. Rare but could happen. 
    # get all ca coords
    ca_coords = pdb_dict['ca_coords']
    IDR_length = len(ca_coords)

    # a few quick checks. 
    if IDR_end_to_end_dist > IDR_length*bond_length:
        raise dodoException('IDR length is not sufficient to be specified IDR_end_to_end_dist!')

    # list to hold IDR coorindates
    IDR_coords=[]

    # end coordinate just going to be diestnace from ca_coords[0] Doesn't really matter tbh.
    end_coord = random_coordinate_on_sphere_surface(ca_coords[0], IDR_end_to_end_dist)
    # add first coord in ca coords to IDR coords. 
    IDR_coords.append(ca_coords[0])

    # Dict of final distances. Basically this lets us specify the sphere radius
    # for the final steps to build the IDR. I determined this more or less empirically
    # and found that without doing this, building the IDR is only successful ~5% of the time.
    final_distances={10: 13, 9: 11.8, 8: 11.6, 7: 11.2, 6: 10.6, 5: 9.8, 4: 9, 3: 8, 2: 6, 1: 3.8}
    if IDR_length > 10:
        radius_distances=sorted(np.linspace(14, IDR_end_to_end_dist, IDR_length-10),reverse=True)
        # finish radius distances
        for i in range(10,0,-1):
            radius_distances.append(final_distances[i])
    else:
        radius_distances=[]
        for i in range(IDR_length,0,-1):
            radius_distances.append(final_distances[i])    

    # dict for distance per residue.
    dist_per_residue = {}
    for i in range(0, IDR_length):
        dist_per_residue[i]=radius_distances[i]
    
    # now can start from start coord and go from there with ever decreasing sphere size.
    for cur_IDR_res in range(1, IDR_length):
        current_radius_dist = dist_per_residue[cur_IDR_res]
        current_start = IDR_coords[-1]
        current_attempt=0
        # keep track of we succeeded. 
        success=False
        for attempt in range(0, attempts_all_coords):
            current_attempt+=1
            new_IDR_coord=random_walk_step(current_start, bond_length)
            if is_min_distance_satisfied(IDR_coords, new_IDR_coord, clash_dist):
                if is_point_inside_sphere(new_IDR_coord, end_coord, current_radius_dist):
                    IDR_coords.append(new_IDR_coord)
                    success=True
                    break
        if success==False:
            raise dodoException(f'Unable to build IDR residue {cur_IDR_res} after {current_attempt} attempts.')

    # final check for bond lengths.
    if bond_length_check(IDR_coords, min_dist=min_bond_dist, 
        max_dist=max_bond_dist)==False:
        raise dodoException('Error! bond_length_check function detected a bond length outside of the min and max bond length range!')

    # verify length is good
    if len(IDR_coords) != len(ca_coords):
        raise dodoException('Error! Length of final coordinates not equal to length of input cooridnates')

    # overwrite final coordinates in pdb dict
    pdb_dict['ca_coords']=IDR_coords
    return pdb_dict
    

def build_terminal_IDR_to_distance(IDR_end_to_end_dist, pdb_dict, regions_dict,
    relative_IDR_location, bond_length=3.8, clash_dist=3.4, attempts_start_coord=10000, 
    attempts_all_coords=30000, min_bond_dist = 2.8, max_bond_dist=4.4, silent=False, debugging=False):
    '''
    function to build an IDR with a specific number of residues
    to be a specific length. Basically uses random walk. Will build
    the IDR from start coordinate to end coordinate for e2e.
    '''
    # make sure relative IDR location is specified.
    if relative_IDR_location not in ['N', 'C']:
        raise dodoException('Please specify relative IDR location as N or C.')

    # get all ca coords
    all_ca_coords = pdb_dict['ca_coords']

    # need to get some info on FD next to the IDR if there is one. 
    all_regions = list(regions_dict.keys())

    if relative_IDR_location=='N':
        IDR_region_index=regions_dict[all_regions[0]]
        end_coordinate = all_ca_coords[IDR_region_index[1]]
        non_IDR_coordinates=all_ca_coords[IDR_region_index[1]:]
        if len(all_regions)>1:
            closest_fd_index=regions_dict[all_regions[1]]
        else:
            raise dodo_exception('Error! No folded domain next to IDR! Use the build_IDR_specific_dist_no_FDs function!')
    else:
        IDR_region_index=regions_dict[all_regions[-1]]
        end_coordinate = all_ca_coords[IDR_region_index[0]-1]
        non_IDR_coordinates=all_ca_coords[:IDR_region_index[0]]
        if len(all_regions)>1:
            closest_fd_index=regions_dict[all_regions[-1]]
        else:
            raise dodo_exception('Error! No folded domain next to IDR! Use the build_IDR_specific_dist_no_FDs function!')



    # get the length of the IDR
    if relative_IDR_location=='N':
        IDR_sequence = pdb_dict['sequence'][IDR_region_index[0]:IDR_region_index[1]]
    else:
        IDR_sequence = pdb_dict['sequence'][IDR_region_index[0]:IDR_region_index[1]+1]
    IDR_length = len(IDR_sequence)

    # a few quick checks. 
    if IDR_end_to_end_dist > IDR_length*bond_length:
        raise dodoException('IDR length is not sufficient to be specified IDR_end_to_end_dist!')

    if len(end_coordinate) !=3:
        raise dodoException('End coordinate must be a 3D coordinate.')

    # now we want to build an IDR that starts / ends in the direction that is on a line
    # from the center of the FD to the end coordinate. The reason for this is that it will
    # avoid weird clashes with the FD.
    # get the center of the FD

    # loops have nested lists, so need to fix that here.     
    if type(closest_fd_index[0])==list:
        closest_fd_index=closest_fd_index[0]

    FD_center = calculate_center_coordinate(all_ca_coords[closest_fd_index[0]:closest_fd_index[1]])
    distance_FD_center_to_start = get_res_dist(FD_center, end_coordinate)
    start_coordinate = generate_coordinate_on_line(FD_center, end_coordinate, IDR_end_to_end_dist+distance_FD_center_to_start)

    # list to hold IDR coorindates
    IDR_coords=[]
    check_these_coordinates=[]
    check_these_coordinates.extend(non_IDR_coordinates)

    # see if there are IDR coordinates that will be reubuilt later anyways that we can ignore
    # this can only really work if we are building the N-terminal IDr. 
    if relative_IDR_location=='N':
        IDRs=[region for region in all_regions if 'idr' in region]
        for idr in IDRs:
            if idr != all_regions[0]:
                idr_coords = all_ca_coords[regions_dict[idr][0]:regions_dict[idr][1]]
                for coord in idr_coords:
                    check_these_coordinates.remove(coord)

    if is_min_distance_satisfied(non_IDR_coordinates, start_coordinate, clash_dist):
        IDR_coords.append(start_coordinate)
        check_these_coordinates.append(start_coordinate)
    else:
        if silent==False:
            print('Falling back to random IDR coordinate start point generation due to clash!')
        # now start trying to get a starting coordinate. 
        cur_attempt=0
        for attempt in range(0, attempts_start_coord):
            cur_attempt+=1
            current_start = random_coordinate_on_sphere_surface(end_coordinate, IDR_end_to_end_dist)
            if is_min_distance_satisfied(non_IDR_coordinates, current_start, clash_dist):
                IDR_coords.append(current_start)
                check_these_coordinates.append(current_start)
                break  

    # Dict of final distances. Basically this lets us specify the sphere radius
    # for the final steps to build the IDR. I determined this more or less empirically
    # and found that without doing this, building the IDR is only successful ~5% of the time.
    final_distances={10: 13, 9: 11.8, 8: 11.6, 7: 11.2, 6: 10.6, 5: 9.8, 4: 9, 3: 8, 2: 6, 1: 3.8}
    if IDR_length > 10:
        radius_distances=sorted(np.linspace(14, IDR_end_to_end_dist, IDR_length-10),reverse=True)
        # finish radius distances
        for i in range(10,0,-1):
            radius_distances.append(final_distances[i])
    else:
        radius_distances=[]
        for i in range(IDR_length,0,-1):
            radius_distances.append(final_distances[i])    

    # dict for distance per residue.
    dist_per_residue = {}
    for i in range(0, IDR_length):
        dist_per_residue[i]=radius_distances[i]
    
    # now can start from start coord and go from there with ever decreasing sphere size.
    for cur_IDR_res in range(1, IDR_length):
        current_radius_dist = dist_per_residue[cur_IDR_res]
        current_start = IDR_coords[-1]
        current_attempt=0
        # keep track of we succeeded. 
        success=False
        for attempt in range(0, attempts_all_coords):
            current_attempt+=1
            new_IDR_coord=random_walk_step(current_start, bond_length)
            if is_min_distance_satisfied(check_these_coordinates, new_IDR_coord, clash_dist):
                if is_point_inside_sphere(new_IDR_coord, end_coordinate, current_radius_dist):
                    IDR_coords.append(new_IDR_coord)
                    check_these_coordinates.append(new_IDR_coord)
                    success=True
                    break
        if success==False:
            raise dodoException(f'Unable to build IDR residue {cur_IDR_res} after {current_attempt} attempts.')

    # rebuild all coordinates.
    final_coordinates=[]
    if relative_IDR_location=='N':
        final_coordinates=IDR_coords
        final_coordinates.extend(non_IDR_coordinates)
    else:
        final_coordinates=non_IDR_coordinates
        rev_IDR_coords=[]
        for i in range(len(IDR_coords)-1,-1,-1):
            rev_IDR_coords.append(IDR_coords[i])
        final_coordinates.extend(rev_IDR_coords)

    # verify length is good
    if len(final_coordinates) != len(all_ca_coords):
        raise dodoException('Error! Length of final coordinates not equal to length of input cooridnates')

    # final check for bond lengths.
    if bond_length_check(final_coordinates, min_dist=min_bond_dist, 
        max_dist=max_bond_dist)==False:
        raise dodoException('Error! bond_length_check function detected a bond length outside of the min and max bond length range!')

    # overwrite final coordinates in pdb dict
    pdb_dict['ca_coords']=final_coordinates
    return pdb_dict
    

def build_fd_loops(pdb_dict, regions_dict, loop_name, bond_length=3.8, clash_dist=3.4,
    attempts_all_coords=30000, min_bond_dist = 2.8, max_bond_dist=4.4, silent=False):
    # function for adding in loops. 

    # get ca coords, loop coords, and fd coords
    ca_coords = pdb_dict['ca_coords']
    loop_coords = regions_dict[loop_name][1]
    fd_coords = regions_dict[loop_name][0]

    # Dict of final distances. Basically this lets us specify the sphere radius
    # for the final steps to build the IDR. I determined this more or less empirically
    # and found that without doing this, building the IDR is only successful ~5% of the time.
    final_distances={10: 21.8, 9: 19.8, 8: 17.8, 7: 15.8, 6: 13.8, 5: 11.8, 4: 9.8, 3: 7.8, 2: 5.8, 1: 3.8}

    # figure out which coords to avoid hitting.
    avoid_these_coords=[]
    for region_name in regions_dict:
        region=regions_dict[region_name]
        if type(region[0])==list:
            all_coords_region = [aa for aa in range(region[0][0],region[0][1])]
            in_fd_loops = region[1]
            for cur_l in in_fd_loops:
                for aa in range(cur_l[0], cur_l[1]):
                    all_coords_region.remove(aa)
            for aa in all_coords_region:
                avoid_these_coords.append(ca_coords[aa])
        else:
            for aa in range(region[0], region[1]):
                avoid_these_coords.append(ca_coords[aa])

    # by residue ca  locs that we can overwrite as we build the loop
    ca_locs_by_res = {}
    for i in range(0, len(ca_coords)):
        ca_locs_by_res[i]=ca_coords[i]
    
    # for each loop
    for cur_loop in loop_coords:
        all_loop_indices = [aa for aa in range(cur_loop[0], cur_loop[1])]
        len_cur_loop = len(all_loop_indices)
        location_of_radius_center = ca_locs_by_res[all_loop_indices[-1]+1]
        max_dist = (len_cur_loop*bond_length)+bond_length
        min_dist = get_res_dist(ca_coords[cur_loop[0]-1], location_of_radius_center)
        if max_dist < min_dist:
            raise dodoException('Number of residues for Loop is insufficient to reach start of loop to end based on bond length')
        
        # get the distances
        if len_cur_loop>10:
            radius_distances=sorted(np.linspace(22, (22+(bond_length*(len_cur_loop-10))), len_cur_loop-10),reverse=True)
            # finish radius len_cur_loop
            for i in range(10,0,-1):
                radius_distances.append(final_distances[i])
        else:
            radius_distances=[]
            for i in range(len(idr_seq),0,-1):
                radius_distances.append(final_distances[i])

        per_res_distances={}
        for radius_ind, loop_res_ind in enumerate(all_loop_indices):
            per_res_distances[loop_res_ind]=radius_distances[radius_ind]
        
        # start to build loop
        cur_attempt=0
        for cur_loop_res in all_loop_indices:
            # track success
            success=False
            cur_dist_from_center = per_res_distances[cur_loop_res]
            start_step_loc = ca_locs_by_res[cur_loop_res-1]
            while cur_attempt < attempts_all_coords:
                cur_attempt+=1
                new_loop_coord=random_walk_step(start_step_loc, bond_length)
                if is_min_distance_satisfied(avoid_these_coords, new_loop_coord, clash_dist):
                    if is_point_inside_sphere(new_loop_coord, location_of_radius_center, cur_dist_from_center):
                        ca_locs_by_res[cur_loop_res]=new_loop_coord
                        avoid_these_coords.append(new_loop_coord)
                        success=True
                        break
            if success==False:
                raise dodoException(f'Unable to build loop residue {cur_loop_res} after {cur_attempt} attempts.')

    # rebuild all coordinates.
    rebuilt_coords = list(ca_locs_by_res.values())

    # check bond lengths
    if bond_length_check(rebuilt_coords, min_dist=min_bond_dist, 
        max_dist=max_bond_dist)==False:
        raise dodoException('Error! bond_length_check function detected a bond length outside of the min and max bond length range!')

    if len(rebuilt_coords) != len(ca_coords):
        raise dodoException('Error! Length of final coordinates not equal to length of input cooridnates')

    # overwrite final coordinates in pdb dict
    pdb_dict['ca_coords']=rebuilt_coords
    return pdb_dict


def build_idr_between_fds(IDR_end_to_end_dist, pdb_dict, regions_dict,
    idr_name, bond_length=3.8, clash_dist=3.4, attempts_all_coords=30000, 
    min_bond_dist = 2.8, max_bond_dist=4.4, silent=False):
    
    # function for building IDRs from FDs. 


    # get some info. 
    ca_coords = pdb_dict['ca_coords']
    idr_indices = regions_dict[idr_name]
    all_idr_indices=[aa for aa in range(idr_indices[0], idr_indices[1])]
    region_names = list(regions_dict.keys())
    idr_index = region_names.index(idr_name)
    first_fd=regions_dict[region_names[idr_index-1]]
    second_fd=regions_dict[region_names[idr_index+1]]
    # check for nested value due to loop.
    if type(first_fd[0])==list:
        first_fd = first_fd[0]
    first_fd_coords=ca_coords[first_fd[0]:first_fd[1]]
    if type(second_fd[0])==list:
        second_fd=second_fd[0]
    second_fd_coords = ca_coords[second_fd[0]:second_fd[1]]
    begin_second_fd_index = second_fd[0]

    # second fd is moved relative to the first so we can get e2e, so
    # we need to move everything from the first res of the second fd on. 
    target_move_region = ca_coords[begin_second_fd_index:]

    # get idr lenght
    idr_length = len([aa for aa in range(idr_indices[0], idr_indices[1])])

    # figure out if it's even theoretically possible to build the IDR.
    if IDR_end_to_end_dist > (idr_length*bond_length): 
        raise dodoException('Cannot build IDR. Distance between FDs is too large.')

    # first thing we need to do is move the second fd relative to the first. 
    # get rough center of the first fd. 
    center_fd1_coord = find_closest_coordinate(first_fd_coords, calculate_center_coordinate(first_fd_coords))
    # get rough center of the second fd. 
    center_fd2_coord = find_closest_coordinate(second_fd_coords, calculate_center_coordinate(second_fd_coords))
    # now make a line along the center coords a specific distance away and get point
    start_coord_second_region = generate_coordinate_on_line(center_fd1_coord, center_fd2_coord, IDR_end_to_end_dist)

    # move the second region
    moved_second_region = translate_coordinates(target_move_region, target_move_region[0], start_coord_second_region)
    # get the radius start point
    radius_start_point=moved_second_region[0]

    # now need to make a new list of the coords by residue.
    ca_coords_after_move = []
    ca_coords_after_move.extend(ca_coords[:begin_second_fd_index])
    ca_coords_after_move.extend(moved_second_region)
    ca_coords_by_res={}
    for i in range(0, len(ca_coords_after_move)):
        ca_coords_by_res[i]=ca_coords_after_move[i]
    
    # list for avoiding coords.
    avoid_coords=[]
    for cur_coord in ca_coords_by_res:
        if cur_coord not in all_idr_indices:
            avoid_coords.append(ca_coords_by_res[cur_coord])

    # Dict of final distances. Basically this lets us specify the sphere radius
    # for the final steps to build the IDR. I determined this more or less empirically
    # and found that without doing this, building the IDR is only successful ~5% of the time.
    final_distances={10: 13, 9: 11.8, 8: 11.6, 7: 11.2, 6: 10.6, 5: 9.8, 4: 9, 3: 8, 2: 6, 1: 3.8}
    if idr_length>10:
        radius_distances=sorted(np.linspace(12, (bond_length*4)+IDR_end_to_end_dist, idr_length-10),reverse=True)
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
    
    
    # Now need to generate locations for each amino acid of the IDR.
    for res in all_idr_indices:
        cur_radius_dist = radius_distances_by_res[res]
        start_res = ca_coords_by_res[res-1]
        cur_attempt=0
        success=False
        while cur_attempt < attempts_all_coords:
            cur_attempt+=1
            new_IDR_coord=random_walk_step(start_res, bond_length)
            if is_min_distance_satisfied(avoid_coords, new_IDR_coord, clash_dist):
                if is_point_inside_sphere(new_IDR_coord, radius_start_point, cur_radius_dist):
                    ca_coords_by_res[res]=new_IDR_coord
                    avoid_coords.append(new_IDR_coord)
                    success=True
                    break
        if success==False:
            raise dodoException(f'Unable to build IDR residue {res} after {cur_attempt} attempts.')
    
    # now some checks.
    rebuilt_coords = list(ca_coords_by_res.values())
    if len(rebuilt_coords)!= len(ca_coords):
        raise dodoException('Error! Length of rebuilt coordinates not equal to length of input cooridnates')

    # check bond lengths
    if bond_length_check(rebuilt_coords, min_dist=min_bond_dist, 
        max_dist=max_bond_dist)==False:
        raise dodoException('Error! bond_length_check function detected a bond length outside of the min and max bond length range!')

    # overwrite coords
    pdb_dict['ca_coords']=rebuilt_coords
    return pdb_dict

def build_structure(pdb_dict, region_dict, mode, attempts_per_region=50,
    attempts_per_coord=5000, bond_length=3.8, clash_dist=3.4, min_bond_dist=2.8, 
    max_bond_dist=4.4, silent=False, verbose=False, debugging=False):

    # final function. 
    run_prediction=False
    if mode =='super_compact':
        length_multiplier = 0.3
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

    idrs = [region for region in region_dict if 'idr' in region]
    fd_loops = [region for region in region_dict if 'loop' in region]
    fds = [region for region in region_dict if 'folded' in region]

    if verbose==True:
        print(f'{len(idrs)} IDRs detected, {len(fd_loops)} FDs with loops detected, {len(fds)} folded domaind detected.')


    # if only IDRs, just do that. 
    if len(fd_loops)==0 and len(fds)==0 and len(idrs) != 0:
        if len(idrs)!=1:
            raise dodoException('Error in /structure_constructor.build_structure() function. Assumed only 1 IDR of there are no FDs or FDs with loops.')
        # if just an idr, shouldn't have multiple...
        idr_seq = pdb_dict['sequence']
        if run_prediction==True:
            e2e = predict_e2e(idr_seq)
        else:
            e2e = len(idr_seq)*length_multiplier
        return build_IDR_specific_dist_no_FDs(e2e, pdb_dict)
    
    else:
        if len(fds)!=0 and len(fd_loops)==0 and len(idrs)==0:
            raise dodoException('Only Fd was input, need an IDR or loop to build.')
        
        # now can get to work. 
        # first build loops
        if len(fd_loops)!= 0:
            for loop in fd_loops:
                success=False
                cur_attempt=0
                while cur_attempt<attempts_per_region:
                    cur_attempt+=1
                    if verbose==True:
                        print(f'On loop {loop}, attempt number {cur_attempt}')
                    try:
                        pdb_dict = build_fd_loops(pdb_dict, region_dict, loop, bond_length=bond_length,
                            clash_dist=clash_dist, attempts_all_coords=attempts_per_coord, 
                            min_bond_dist=min_bond_dist, max_bond_dist=max_bond_dist, silent=silent)
                        success=True
                        break
                    except dodoException:
                        pass
                if success==False:
                    if debugging==True:
                        print(region_dict)
                        print(loop)
                    raise dodoException(f'Unable to build loop {loop} after {cur_attempt} attempts.')

        # now for the IDRs. 
        if len(idrs)!=0:
            terminal_IDRs={}
            for IDR in idrs:
                if IDR == list(region_dict)[0]:
                    terminal_IDRs['N']=IDR
                elif IDR == list(region_dict)[-1]:
                    terminal_IDRs['C']=IDR
                else:
                    # get IDR length
                    IDR_seq = pdb_dict['sequence'][region_dict[IDR][0]:region_dict[IDR][1]]
                    if run_prediction==True:
                        e2e = predict_e2e(IDR_seq)
                    else:
                        e2e=len(IDR_seq)*length_multiplier
                    success=False
                    cur_attempt=0
                    while cur_attempt<attempts_per_region:
                        cur_attempt+=1
                        if verbose==True:
                            print(f'On IDR {IDR}, attempt number {cur_attempt}')
                        try:
                            pdb_dict = build_idr_between_fds(e2e, pdb_dict, region_dict, IDR, 
                                bond_length=bond_length, clash_dist=clash_dist, 
                                attempts_all_coords=attempts_per_coord, 
                                min_bond_dist=min_bond_dist, max_bond_dist=max_bond_dist, silent=silent)
                            success=True
                            break
                        except dodoException:
                            pass
                    if success==False:
                        if debugging==True:
                            print(region_dict)
                            print(IDR_seq)
                        raise dodoException(f'Unable to build IDR {IDR} after {cur_attempt} attempts.')
            # now for the terminal IDRs
            if 'N' in list(terminal_IDRs):
                # get IDR length
                IDR = terminal_IDRs['N']
                IDR_seq = pdb_dict['sequence'][region_dict[IDR][0]:region_dict[IDR][1]]
                if run_prediction==True:
                    e2e = predict_e2e(IDR_seq)
                else:
                    e2e=len(IDR_seq)*length_multiplier
                success=False
                cur_attempt=0
                while cur_attempt<attempts_per_region:
                    cur_attempt+=1
                    if verbose==True:
                        print(f'On N-terminal IDR, attempt number {cur_attempt}')
                    try:
                        pdb_dict = build_terminal_IDR_to_distance(e2e, pdb_dict, region_dict, 'N', 
                            bond_length=bond_length, clash_dist=clash_dist, 
                            attempts_start_coord=20000, attempts_all_coords=attempts_per_coord, 
                            min_bond_dist=min_bond_dist, max_bond_dist=max_bond_dist, silent=silent)
                        success=True
                        break
                    except dodoException:
                        pass
                if success==False:
                    if debugging==True:
                        print(region_dict)
                        print(IDR_seq)                    
                    raise dodoException(f'Unable to build IDR {IDR} after {cur_attempt} attempts.')                
            if 'C' in list(terminal_IDRs):
                # get IDR length
                IDR = terminal_IDRs['C']
                IDR_seq = pdb_dict['sequence'][region_dict[IDR][0]:region_dict[IDR][1]]
                if run_prediction==True:
                    e2e = predict_e2e(IDR_seq)
                else:
                    e2e=len(IDR_seq)*length_multiplier
                success=False
                cur_attempt=0
                while cur_attempt<attempts_per_region:
                    cur_attempt+=1
                    if verbose==True:
                        print(f'On N-terminal IDR, attempt number {cur_attempt}')
                    try:
                        pdb_dict = build_terminal_IDR_to_distance(e2e, pdb_dict, region_dict, 'C', 
                            bond_length=bond_length, clash_dist=clash_dist, 
                            attempts_start_coord=20000, attempts_all_coords=attempts_per_coord, 
                            min_bond_dist=min_bond_dist, max_bond_dist=max_bond_dist, silent=silent)
                        success=True
                        break
                    except dodoException:
                        pass
                if success==False:
                    if debugging==True:
                        print(region_dict)
                        print(IDR_seq)                         
                    raise dodoException(f'Unable to build IDR {IDR} after {cur_attempt} attempts.') 
    return pdb_dict
