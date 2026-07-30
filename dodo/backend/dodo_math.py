'''
holds common math and math-related functions.
'''
import numpy as np

def calculate_distance(coord1, coord2):
    '''
    This function calculates the Euclidean distance between two points.
    
    Parameters
    ----------
    coord1 : np.array, shape=(3,)
        coordinates of the first point
    coord2 : np.array, shape=(3,)
        coordinates of the second point
    
    Returns
    -------
    float
        Euclidean distance between the two points
    '''
    return np.linalg.norm(coord1 - coord2)

def calculate_distance_vectorized(coord_list1, coord_list2):
    '''
    This function calculates the Euclidean distance between two lists of points.
    
    Parameters
    ----------
    coord_list1 : np.array, shape=(N, 3)
        coordinates of the first list of points
    coord_list2 : np.array, shape=(M, 3)
        coordinates of the second list of points
    
    Returns
    -------
    np.array, shape=(N, M)
        Euclidean distance between each pair of points
    '''
    return np.linalg.norm(coord_list1[:, np.newaxis, :] - coord_list2[np.newaxis, :, :], axis=2)



def find_furthest_coordinate(coord_list1, coord_list2):
    '''
    Finds the coordinate in the first list that is furthest from all coordinates in the second list
    using vectorized operations for efficiency.
    
    Parameters
    ----------
    coord_list1 : array-like, shape=(N, 3)
        First list of coordinates to check
    coord_list2 : array-like, shape=(M, 3)
        Second list of coordinates to compare against
    
    Returns
    -------
    furthest_coord : numpy.ndarray, shape=(3,)
        The coordinate from coord_list1 that has the maximum minimum distance to any point in coord_list2
    max_min_distance : float
        The minimum distance from the returned coordinate to any point in coord_list2
    
    Raises
    ------
    ValueError
        If inputs have incorrect shapes or are empty
    '''
    try:
        points1 = np.asarray(coord_list1, dtype=float)
        points2 = np.asarray(coord_list2, dtype=float)
    except (ValueError, TypeError):
        raise ValueError("Coordinates must be convertible to float arrays")

    # Validate input shapes
    if points1.ndim != 2 or points1.shape[1] != 3:
        raise ValueError("coord_list1 must be an Nx3 array")
    if points2.ndim != 2 or points2.shape[1] != 3:
        raise ValueError("coord_list2 must be an Mx3 array")
    if len(points1) == 0 or len(points2) == 0:
        raise ValueError("Coordinate lists cannot be empty")

    # Vectorized distance calculation
    # Reshape arrays for broadcasting: (N, 1, 3) - (1, M, 3) -> (N, M, 3)
    diff = points1[:, np.newaxis, :] - points2[np.newaxis, :, :]
    
    # Calculate distances efficiently using einsum
    # equivalent to np.sqrt(np.sum(diff * diff, axis=2))
    distances = np.sqrt(np.einsum('ijk,ijk->ij', diff, diff))
    
    # Find minimum distance for each point in points1
    min_distances = np.min(distances, axis=1)
    
    # Find the point with the maximum minimum distance
    max_min_idx = np.argmax(min_distances)
    
    return points1[max_min_idx]




def find_transform(coord1, coord2):
    '''
    This function finds the translation and rotation that
    aligns coord1 to coord2.
    '''
    p1 = np.asarray(coord1)
    p2 = np.asarray(coord2)
    
    translation = -p1
    p2_translated = p2 + translation
    
    # Check if points are effectively the same
    norm = np.linalg.norm(p2_translated)
    if np.isclose(norm, 0.0):
        return translation, np.eye(3)
    
    # Calculate rotation matrix to align p2 with x-axis
    v = p2_translated / norm
    x_axis = np.array([1., 0., 0.])
    
    # Optimize rotation calculation
    cos_theta = np.dot(v, x_axis)
    if np.isclose(cos_theta, 1.0):
        return translation, np.eye(3)
    if np.isclose(cos_theta, -1.0):
        return translation, -np.eye(3)
    
    rotation_axis = np.cross(v, x_axis)
    rotation_axis /= np.linalg.norm(rotation_axis)
    
    # Pre-compute values for Rodrigues formula
    sin_theta = np.sqrt(1 - cos_theta**2)
    K = np.array([[0, -rotation_axis[2], rotation_axis[1]],
                  [rotation_axis[2], 0, -rotation_axis[0]],
                  [-rotation_axis[1], rotation_axis[0], 0]])
    
    rotation = np.eye(3) + sin_theta * K + (1 - cos_theta) * (K @ K)
    return translation, rotation

def apply_transform(coord, translation, rotation):
    return rotation @ (np.asarray(coord) + translation)

def reverse_transform(coord, translation, rotation):
    return (rotation.T @ np.asarray(coord)) - translation




def generate_circle_points(origin, center, radius, num_points):
    origin = np.asarray(origin)
    center = np.asarray(center)
    
    normal = center - origin
    normal /= np.linalg.norm(normal)
    
    # Pre-compute angles for all points at once
    angles = np.linspace(0, 2*np.pi, num_points, endpoint=False)
    
    # Optimize orthogonal vector calculation
    v1 = np.array([1., 0., 0.]) if not np.allclose(normal, [1., 0., 0.]) else np.array([0., 1., 0.])
    u = np.cross(normal, v1)
    u /= np.linalg.norm(u)
    v = np.cross(normal, u)
    
    # Vectorized point generation
    cos_theta = np.cos(angles)
    sin_theta = np.sin(angles)
    
    # Generate all points at once using broadcasting
    points = center + radius * (u[:, np.newaxis] * cos_theta + v[:, np.newaxis] * sin_theta).T
    return points

def generate_circle_points_vectorized(origins, centers, radii, num_points):
    origins = np.asarray(origins)
    centers = np.asarray(centers)
    radii = np.asarray(radii)

    normals = centers - origins
    normals /= np.linalg.norm(normals, axis=1, keepdims=True)

    angles = np.linspace(0, 2*np.pi, num_points, endpoint=False)

    v1 = np.array([1., 0., 0.])
    v1 = np.where(np.allclose(normals, v1), np.array([0., 1., 0.]), v1)

    u = np.cross(normals, v1)
    u /= np.linalg.norm(u, axis=1, keepdims=True)
    v = np.cross(normals, u)

    cos_theta = np.cos(angles)
    sin_theta = np.sin(angles)

    points = centers[:, np.newaxis, :] + radii[:, np.newaxis, np.newaxis] * (u[:, np.newaxis, :] * cos_theta + v[:, np.newaxis, :] * sin_theta)
    return points




def generate_random_points_by_distance(origins, distances):
    '''
    Generates random points that are given distances away from the origins.
    
    Parameters
    ----------
    origins : np.array, shape=(N, 3)
        coordinates of the origins
    distances : np.array, shape=(N,)
        distances from the origins

    Returns
    -------
    np.array, shape=(N, 3)
        coordinates of the random points
    '''
    # Generate random points on a unit sphere
    points = np.random.normal(0, 1, (len(origins), 3))
    points /= np.linalg.norm(points, axis=1, keepdims=True)

    # Scale the points to the desired distances
    points *= distances[:, np.newaxis]

    # Translate the points to the origins
    points += origins

    return points




def generate_random_points_by_distance_mult(origins, distances, num_points):
    '''
    Generates multiple random points that are given distances away from the origins.
    
    Parameters
    ----------
    origins : np.array, shape=(N, 3)
        coordinates of the origins
    distances : np.array, shape=(N,)
        distances from the origins
    num_points : int
        number of random points to generate per origin

    Returns
    -------
    np.array, shape=(N, num_points, 3)
        coordinates of the random points
    '''
    N = len(origins)
    points = np.zeros((N, num_points, 3))
    
    for i in range(N):
        # Generate random points on a unit sphere
        rand_points = np.random.normal(0, 1, (num_points, 3))
        rand_points /= np.linalg.norm(rand_points, axis=1, keepdims=True)

        # Scale the points to the desired distances
        rand_points *= distances[i]

        # Translate the points to the origins
        rand_points += origins[i]
        
        points[i] = rand_points

    return points



def find_points_within_sphere(points_of_interest, sphere_center, sphere_radius):
    '''
    Find all points that lie within a sphere of given radius and center.
    
    Parameters
    ----------
    points_of_interest : array-like, shape=(N, 3)
        Array of points to check. Each point should be a 3D coordinate.
    sphere_center : array-like, shape=(3,)
        Center coordinates of the sphere
    sphere_radius : float
        Radius of the sphere. Must be positive.
    
    Returns
    -------
    numpy.ndarray
        Array containing only the points that lie within the sphere
    points_mask : numpy.ndarray
        Boolean mask indicating which points are within the sphere
    
    Raises
    ------
    ValueError
        If inputs have incorrect shapes or if radius is not positive
    '''
    # Convert inputs to numpy arrays and validate
    try:
        points = np.asarray(points_of_interest, dtype=float)
        center = np.asarray(sphere_center, dtype=float)
        radius = float(sphere_radius)
    except (ValueError, TypeError):
        raise ValueError("All inputs must be convertible to float types")

    # Check shapes
    if points.ndim != 2 or points.shape[1] != 3:
        raise ValueError("points_of_interest must be an Nx3 array")
    if center.shape != (3,):
        raise ValueError("sphere_center must be a 3D point")
    if radius <= 0:
        raise ValueError("sphere_radius must be positive")

    # Vectorized distance calculation
    # Using einsum for better memory efficiency with large arrays
    # equivalent to: distances = np.sqrt(np.sum((points - center)**2, axis=1))
    distances = np.sqrt(np.einsum('ij,ij->i', points - center, points - center))
    
    # Return both the points and the mask
    return points[distances <= radius]



def find_point_closest_to_sphere_surface(points_of_interest, 
                                         sphere_center, 
                                         sphere_radius):
    '''
    This function finds the point closest to the surface of a sphere.
    
    Parameters
    ----------
    points_of_interest : np.array, shape=(N, 3)
        points we want to check if they are within the sphere
    sphere_center : np.array, shape=(3,)
        center of the sphere
    sphere_radius : float
        radius of the sphere
    
    Returns
    -------
    np.array, shape=(3,)
        point closest to the sphere surface
    '''
    # calculate distances
    distances = np.linalg.norm(points_of_interest - sphere_center, axis=1)

    # find the point closest to the sphere surface
    closest_point = points_of_interest[np.argmin(np.abs(distances - sphere_radius))]
    return closest_point


def find_points_closest_to_sphere_surface(points_of_interest, 
                                         sphere_center, 
                                         sphere_radius, 
                                         allowed_error=0.2,
                                         return_best_point=True):
    '''
    This function finds the point closest to the surface of a sphere.
    
    Parameters
    ----------
    points_of_interest : np.array, shape=(N, 3)
        points we want to check if they are within the sphere
    sphere_center : np.array, shape=(3,)
        center of the sphere
    sphere_radius : float
        radius of the sphere
    allowed_error : float
        allowed error in the distance from the sphere surface
        this lets us get points close to the surface but within some error
    return_best_point : bool
        if True, return the best point. if False, return all points within the allowed error
    
    Returns
    -------
    np.array, shape=(3,)
        point closest to the sphere surface
    '''
    # calculate distances
    distances = np.linalg.norm(points_of_interest - sphere_center, axis=1)

    # find points within allowed error. 
    # this lets us get points close to the surface but within some error
    points_within_error = points_of_interest[np.abs(distances - sphere_radius) <= allowed_error]
    if len(points_within_error) == 0:
        if return_best_point:
            return [points_of_interest[np.argmin(np.abs(distances - sphere_radius))]]
        else:
            raise ValueError('No points found within allowed error')
    return points_within_error


def get_random_coordinates_on_sphere_surface(starting_coordinate, radius, num_coordinates):
    '''
    Returns a specified number of random coordinates 
    that are uniformly distributed on the surface of a sphere
    with a specified radius and center.
    
    Parameters
    ----------
    starting_coordinate : np.array
        The center of the sphere
    radius : float
        The radius of the sphere
    num_coordinates : int
        The number of random coordinates to return
    
    Returns
    -------
    coordinates : np.array
        An array of coordinates on the surface of the sphere
    '''
    starting_coordinate = np.array(starting_coordinate)
    coordinates = []
    
    for _ in range(num_coordinates):
        # Generate uniform random points using Gaussian distribution
        # This ensures uniform distribution on sphere surface
        x, y, z = np.random.normal(0, 1, 3)
        # Normalize to get point on unit sphere
        vec = np.array([x, y, z])
        vec = vec / np.linalg.norm(vec)
        # Scale by radius and add center coordinates
        point = starting_coordinate + radius * vec
        coordinates.append(point)
    return np.array(coordinates)


def find_points_not_clashing(potential_coords, coords_to_check, clash_distance=3):
    '''
    This function finds the points in potential_coords that are not clashing with any points in coords_to_check.
    
    Parameters
    ----------
    potential_coords : np.array, shape=(N, P, 3)
        Points to check for clashes. N is number of sets, P is points per set.
    coords_to_check : np.array, shape=(M, 3)
        Points to check against.
    clash_distance : float
        Distance at which two points are considered clashing.
    
    Returns
    -------
    np.array, shape=(K, 3)
        Points that are not clashing.
    '''
    # Reshape arrays for broadcasting
    potential_coords = np.asarray(potential_coords)
    coords_to_check = np.asarray(coords_to_check)
    
    # Calculate distances between all points
    # Reshape arrays to allow broadcasting
    distances = np.linalg.norm(
        potential_coords[..., np.newaxis, :] - coords_to_check, 
        axis=-1
    )
    
    # Find points that don't clash with any point in coords_to_check
    non_clashing_mask = np.all(distances > clash_distance, axis=-1)
    non_clashing_points = potential_coords[non_clashing_mask]
    
    return non_clashing_points



def find_points_by_angle(points_to_check, coord_1, coord_2, angle):
    '''
    filter points by being within some angle
    formed by the angle formed by coord_1, coord_2 and points_to_check
    
    Parameters
    ----------
    points_to_check : np.array, shape=(N, 3)
        points we want to check if they are within the angle
    coord_1 : np.array, shape=(3,)
        first coordinate
    coord_2 : np.array, shape=(3,)
        second coordinate
    angle : float
        allowed angle in degrees

    Returns
    -------
    np.array, shape=(K, 3)
        points that are within the allowed angle
    '''
    # Convert inputs to numpy arrays
    points = np.asarray(points_to_check)
    c1 = np.asarray(coord_1)
    c2 = np.asarray(coord_2)
    
    # Convert allowed angle to radians
    allowed_rad = np.radians(angle)
    
    # Calculate vectors from coord_1 to points and coord_2
    v1 = points - c1
    v2 = c2 - c1
    
    # Normalize vectors
    v1_norm = np.linalg.norm(v1, axis=1, keepdims=True)
    v2_norm = np.linalg.norm(v2)
    
    # Avoid division by zero
    mask = v1_norm[:, 0] > 0
    v1_normalized = np.zeros_like(v1)
    v1_normalized[mask] = v1[mask] / v1_norm[mask]
    v2_normalized = v2 / v2_norm
    
    # Calculate angles using dot product
    # Use einsum for efficient dot product calculation
    cos_angles = np.einsum('ij,j->i', v1_normalized, v2_normalized)
    
    # Clip values to avoid numerical errors
    cos_angles = np.clip(cos_angles, -1.0, 1.0)
    angles = np.arccos(cos_angles)
    
    # Return points within allowed angle
    return points[angles <= allowed_rad]


def get_starting_points(num_starting_points, starting_coordinate=[0,0,0], distance=3.8):
    '''
    get an arbitrary number of starting poinst
    distance away from starting coordinate.
    
    Parameters
    ----------
    num_starting_points : int
        number of starting points to generate
    starting_coordinate : np.array, shape=(3,)
        starting coordinate
    
    Returns
    -------
    np.array, shape=(num_starting_points, 2, 3)
        array of starting points
    '''
    # Convert starting_coordinate to numpy array
    starting_coordinate = np.array(starting_coordinate)
    
    # Generate coordinates on sphere surface
    coords = get_random_coordinates_on_sphere_surface(starting_coordinate, distance, num_starting_points)
    
    # Create output array with shape (num_starting_points, 2, 3)
    # Each point pair consists of the starting coordinate and a random point
    result = np.zeros((num_starting_points, 2, 3))
    result[:, 0] = starting_coordinate
    result[:, 1] = coords
    
    return np.array(result)

def get_ending_points(start_points, distances, 
                      clash_points=None, clash_dist=3,
                      num_points_to_check=1000):
    '''
    Get ending points based on start points, distance, and any points to 
    check clashes with. 
    This is used to generate IDRs.

    Parameters
    ----------
    start_points : np.array, shape=(N, 3)
        array of starting points
    distances : np.array, shape=(N,)
        array of distances that generated ending points
        will be from the starting points
    clash_points : np.array, shape=(M, 3)
        array of points to check for clashes
    clash_dist : float
        the distance to use to consider something as clashing with something else. 
    num_points_to_check : int
        number of points to generate per end point that will be checked for clashes.
    
    Returns
    -------
    ending_points : np.array, shape=(N, 3)
        array of ending points
    '''
    # Extract origins from start_points
    origins = start_points
    N = len(start_points)
    
    # Initialize an array to hold the final ending points
    ending_points = np.zeros((N, 3))

    # Generate points for each origin separately
    for i in range(N):
        # Generate multiple random points for this origin
        origin_repeated = np.tile(origins[i], (num_points_to_check, 1))
        potential_coords = generate_random_points_by_distance(
            origin_repeated, 
            np.full(num_points_to_check, distances[i])
        )

        if clash_points is not None:
            # Reshape potential coords for clash checking
            potential_coords_reshaped = potential_coords[np.newaxis, :, :]
            clash_points_reshaped = clash_points[np.newaxis, :, :]

            # Find non-clashing points
            non_clashing_coords = find_points_not_clashing(
                potential_coords_reshaped,
                clash_points_reshaped,
                clash_dist
            )

            # If there are any non-clashing points, select the first one
            if len(non_clashing_coords) > 0:
                ending_points[i] = non_clashing_coords[0]
            else:
                # If all points clash, use the first potential coordinate
                ending_points[i] = potential_coords[0]
        else:
            # If no clash points are provided, use the first generated point
            ending_points[i] = potential_coords[0]

    return ending_points
