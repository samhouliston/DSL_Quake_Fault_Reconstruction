import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import plotly.graph_objs as go
import random
import math
import argparse
import csv

def random_unit_vector():
    # Generate random spherical coordinates
    theta = random.uniform(0, 2*math.pi)  # Azimuthal angle
    phi = random.uniform(0, math.pi)      # Polar angle

    # Convert spherical coordinates to Cartesian coordinates
    x = math.sin(phi) * math.cos(theta)
    y = math.sin(phi) * math.sin(theta)
    z = math.cos(phi)

    return (x, y, z)


def plane_main_axes_from_normal(normal): 
    # Find two axes of plane defined by the center and its normal
    length_axis = random_unit_vector()
    while np.dot(length_axis, normal) < 0.01:
        length_axis = random_unit_vector()

    length_axis = np.cross(normal, length_axis)
    length_axis = length_axis / np.linalg.norm(length_axis)
    width_axis = np.cross(normal, length_axis)

    return length_axis, width_axis


def sample_vec_close_to(target_vec):
    # Sample random vector orthogoanl to target_vec: 
    rand_vec = random_unit_vector()
    rand_vec_orthogonal = np.cross(rand_vec, target_vec)

    # Sample angle between -45 and 45 degrees
    angle = np.random.uniform(- np.pi / 6, np.pi / 6)
    alpha = np.cos(angle)
    beta = np.sin(angle)

    new_vec = alpha * target_vec + beta * rand_vec_orthogonal 
    new_vec = new_vec / np.linalg.norm(new_vec)

    return new_vec 


def sample_length_width_axes_from_favoured_axes(favoured_axis_1, favoured_axis_2):
    favoured_normal = np.cross(favoured_axis_1, favoured_axis_2)

    # SAMPLE THE LENGTH AXIS
    # Choose length axis close to either favoured_axis_1 or favoured_axis_2
    idx = random.randint(0, 1)
    if idx == 0: 
        length_axis = sample_vec_close_to(favoured_axis_1)
    else: 
        length_axis = sample_vec_close_to(favoured_axis_2)

    # OBTRAIN WIDHT AXIS
    width_axis = np.cross(length_axis, favoured_normal)
    normal_axis = np.cross(length_axis, width_axis)
    
    return length_axis, width_axis, normal_axis 



def display_fault_properties(fault_metadata):
    (center, length, width, normal, length_axis, width_axis) = fault_metadata
    print(f"         center: {center} \n         length: {length}, width: {width} \n         length_axis: {length_axis},   width_axis: {width_axis}")


def sample_new_center(centers_list, favoured_axis_1, favoured_axis_2, domain_scale=1, eps_center=0.6):
    """
    Samples a new center at least eps_center distance from previous centers in centers_list.
    The new center is sampled either preferable along axis_1 or axis_2. 
    """
    third_axis = np.cross(favoured_axis_1, favoured_axis_2)
    length = domain_scale

    # Ensure centers aren't too close
    while True: 
            # Choose to favour one axis or the other: 
            idx = random.randint(1,2)
            if idx == 1:
                alpha = np.random.uniform(-2 *length, 2 *length)
                beta = np.random.uniform(-length/3, length/3)
                gamma = np.random.uniform(-length/3, length/3)
            
            else:
                alpha = np.random.uniform(-length/3, length/3)
                beta = np.random.uniform(-2*length, 2*length)
                gamma = np.random.uniform(-length/3, length/3)
                
            center = alpha * favoured_axis_1 + beta * favoured_axis_2 + gamma * third_axis
            np.random.uniform(-domain_scale, domain_scale, size=3)
            
            if all(np.linalg.norm(center - c) > 2 * eps_center for c in centers_list):  
                break

    return center



def generate_faults_metadata(n_simple_faults, n_bent_faults, n_cross_faults, n_Y_faults, n_parallel_faults, n_ladder_structures, eps_center=0.5, domain_scale=1, VERBOSE=False):
    """
    Generate metadata characterizing faults via a plane defined by a center, length, width, normal, length_axis, and width_axis.
    
    Args:
    - n_simple_faults:      number of simple planar faults to generate
    - n_bent_faults:        number of bent (bi-planar) faults to generate
    - n_cross_faults:       number of bent (bi-planar) faults to generate
    - n_Y_faults:           number of bent (bi-planar) faults to generate
    - n_parallel_faults:    number of bent (bi-planar) faults to generate
    - n_ladder_structures:  number of ladder-like fault structures to generate
    - eps_center:           minimum distance between fault centers
    - domain_scale:         scale=x denotes [-x,x]^3 space. 
    
    Returns:
    - faults_info: List of tuples containing fault metadata (center, length, width, normal, length_axis, width_axis)
    """
    
    # Set fault plane max length and width as proportion of domain
    L_max = 1.5 * domain_scale
    L_min = 1 * domain_scale
    W_max = 0.7 * domain_scale
    W_min= 0.5 * domain_scale
    faults_info = []
    centers_list = []
    # Generate two orthogonal directions to align all faults
    random_unit_vec = random_unit_vector()
    favoured_axis_1, favoured_axis_2 = plane_main_axes_from_normal(random_unit_vec)


    # Generate simple fault data
    if VERBOSE:
        print(f"--------------------------------------------------")
        print(f"    Generating {n_simple_faults} simple faults:")
        print(f"--------------------------------------------------")
    for i_simp in range(n_simple_faults):    
        center = sample_new_center(centers_list, favoured_axis_1, favoured_axis_2, domain_scale=1, eps_center=eps_center)
        centers_list.append(center)
        length = np.random.uniform(L_min, L_max)
        width = np.random.uniform(W_min, W_max)
        #normal = random_unit_vector()
        #length_axis, width_axis = plane_main_axes_from_normal(normal)
        length_axis, width_axis, normal = sample_length_width_axes_from_favoured_axes(favoured_axis_1, favoured_axis_2)
        faults_info.append((center, length, width, normal, length_axis, width_axis))

        if VERBOSE:
            print(f"    Simple fault number: {i_simp+1}")
            display_fault_properties((center, length, width, normal, length_axis, width_axis))


    # Generate bent fault data
    if VERBOSE:
        print(f"\n--------------------------------------------------")
        print(f"    Generating {n_bent_faults} bent faults:")
        print(f"--------------------------------------------------")
    for i in range(n_bent_faults):

        # Generate the main fault data
        main_center = sample_new_center(centers_list, favoured_axis_1, favoured_axis_2, domain_scale=1, eps_center=eps_center)
        centers_list.append(main_center)
        main_length = np.random.uniform(L_min, L_max)
        main_width = np.random.uniform(W_min, W_max)
        #main_normal = random_unit_vector()
        #main_length_axis, main_width_axis = plane_main_axes_from_normal(main_normal)
        main_length_axis, main_width_axis, main_normal = sample_length_width_axes_from_favoured_axes(favoured_axis_1, favoured_axis_2)
        faults_info.append((main_center, main_length, main_width, main_normal, main_length_axis, main_width_axis))

        if VERBOSE:
            print(f"    Bent fault number: {i//2 + 1}")
            print("       First plane:")
            display_fault_properties((main_center, main_length, main_width, main_normal, main_length_axis, main_width_axis))

        # Generate secondary smaller fault: sample plane at an angle
        bent_factor = np.random.uniform(0.6, 0.75)
        bent_length = bent_factor * main_length
        bent_width = bent_factor * main_width 

        # Sample an angle, calculate new axes
        angle = np.random.uniform(-50, 50, size=1)
        alpha = np.cos(np.radians(90 - angle))
        beta = np.sqrt(1-alpha**2)
        bent_length_axis = beta*np.array(main_length_axis)  + alpha*np.array(main_normal)
        bent_length_axis = bent_length_axis / np.linalg.norm(bent_length_axis)
        bent_width_axis = main_width_axis 
        
        bent_center = main_center + 0.5*main_length*main_length_axis + 0.5*bent_length * bent_length_axis
        bent_normal = np.cross(bent_length_axis, bent_width_axis)

        faults_info.append((bent_center, bent_length, bent_width, bent_normal, bent_length_axis, bent_width_axis))

        if VERBOSE:
            print(f"\n       Second plane:")
            print(f"         Angle: {angle}")
            display_fault_properties((bent_center, bent_length, bent_width, bent_normal, bent_length_axis, bent_width_axis))


    # Generate cross fault data
    if VERBOSE:
        print(f"\n--------------------------------------------------")
        print(f"    Generating {n_cross_faults} cross faults:")
        print(f"--------------------------------------------------")
    for k in range(n_cross_faults):

        # Generate main fault data
        cross_center = sample_new_center(centers_list, favoured_axis_1, favoured_axis_2, domain_scale=1, eps_center=eps_center)
        centers_list.append(cross_center)
        main_length = np.random.uniform(L_min, L_max)
        main_width = np.random.uniform(W_min, W_max)
        #main_normal = random_unit_vector()
        #main_length_axis, main_width_axis = plane_main_axes_from_normal(main_normal)
        main_length_axis, main_width_axis, main_normal = sample_length_width_axes_from_favoured_axes(favoured_axis_1, favoured_axis_2)
        faults_info.append((cross_center, main_length, main_width, main_normal, main_length_axis, main_width_axis))

        if VERBOSE:
            print(f"    Cross fault number: {k//2 + 1}")
            print("       First plane:")
            display_fault_properties((cross_center, main_length, main_width, main_normal, main_length_axis, main_width_axis))

       # Generate secondary perpendicular 'crossed' fault
        cross_factor = 1
        crossed_length = cross_factor * main_length
        crossed_width = 0.8*cross_factor * main_width 
        crossed_normal = main_length_axis
        crossed_length_axis = main_normal 
        crossed_width_axis = np.cross(crossed_length_axis, crossed_normal)
    
        faults_info.append((cross_center, crossed_length, crossed_width, crossed_normal, crossed_length_axis, crossed_width_axis))

        if VERBOSE:
            print(f"\n       Crosssed plane:")
            display_fault_properties((cross_center, crossed_length, crossed_width, crossed_normal, crossed_length_axis, crossed_width_axis))


    # Generate Y-shaped fault data
    if VERBOSE:
        print(f"\n--------------------------------------------------")
        print(f"    Generating {n_Y_faults} Y - faults:")
        print(f"--------------------------------------------------")
    for l in range(n_Y_faults):

        # Generate main fault data
        main_center = sample_new_center(centers_list, favoured_axis_1, favoured_axis_2, domain_scale=1, eps_center=eps_center)
        centers_list.append(main_center)
        main_length = 1.4 * np.random.uniform(L_min, L_max)
        main_width = 0.7* np.random.uniform(W_min, W_max)
        #main_normal = random_unit_vector()
        #main_length_axis, main_width_axis = plane_main_axes_from_normal(main_normal)
        main_length_axis, main_width_axis, main_normal = sample_length_width_axes_from_favoured_axes(favoured_axis_1, favoured_axis_2)
        faults_info.append((main_center, main_length, main_width, main_normal, main_length_axis, main_width_axis))

        if VERBOSE:
            print(f"    Y - fault number: {l//2 + 1}")
            print("       First plane:")
            display_fault_properties((main_center, main_length, main_width, main_normal, main_length_axis, main_width_axis))

        # Generate secondary perpendicular 'branched' fault
        branch_length = 0.35*main_length
        branch_width = 0.7*main_width 

        angle = np.random.choice([-1, 1]) * np.random.uniform(20, 35, size=1)
        alpha = np.cos(np.radians(90 - angle))
        beta = np.sqrt(1-alpha**2)
        branch_length_axis = beta*np.array(main_length_axis)  + alpha*np.array(main_normal)
        branch_length_axis = branch_length_axis / np.linalg.norm(branch_length_axis)
        branch_width_axis = main_width_axis    
        branch_normal = np.cross(branch_length_axis, branch_width_axis)    

        branch_center = main_center + np.random.uniform(0, 0.15)*main_length*main_length_axis + 0.5* branch_length * branch_length_axis
   
        
        faults_info.append((branch_center, branch_length, branch_width, branch_normal, branch_length_axis, branch_width_axis))

        if VERBOSE:
            print(f"\n       Branch plane:")
            display_fault_properties((branch_center, branch_length, branch_width, branch_normal, branch_length_axis, branch_width_axis))


    # Generate two parallel fault data
    if VERBOSE:
        print(f"\n--------------------------------------------------")
        print(f"    Generating {n_parallel_faults} parallel faults:")
        print(f"--------------------------------------------------")
    for b in range(n_parallel_faults):

        # Generate main fault data
        main_center = sample_new_center(centers_list, favoured_axis_1, favoured_axis_2, domain_scale=1, eps_center=eps_center)
        centers_list.append(main_center)
        main_length = np.random.uniform(L_min, L_max)
        main_width = np.random.uniform(W_min, W_max)
        #main_normal = random_unit_vector()
        #main_length_axis, main_width_axis = plane_main_axes_from_normal(main_normal)
        main_length_axis, main_width_axis, main_normal = sample_length_width_axes_from_favoured_axes(favoured_axis_1, favoured_axis_2)
        faults_info.append((main_center, main_length, main_width, main_normal, main_length_axis, main_width_axis))

        if VERBOSE:
            print(f"    Parallel pair number: {b//2 + 1}")
            print("       First plane:")
            display_fault_properties((main_center, main_length, main_width, main_normal, main_length_axis, main_width_axis))

        # Generate secondary perpendicular 'crossed' fault
        dist_between_parallels = np.random.uniform(0.5*main_width, 0.7*main_width)
        distance_shifted = np.random.uniform(0.2, 0.7) * main_length
        second_center = main_center + dist_between_parallels * np.array(main_normal) + distance_shifted * main_length_axis
        second_factor = np.random.uniform(0.6, 1.5)
        second_length = second_factor * main_length
        second_width = second_factor * main_width 
        second_length_axis = main_length_axis
        second_width_axis = main_width_axis
        second_normal = np.cross(second_length_axis, second_width_axis)
        faults_info.append((second_center, second_length, second_width, second_normal, second_length_axis, second_width_axis))

        if VERBOSE:
            print(f"\n       Second plane:")
            display_fault_properties((second_center, second_length, second_width, second_normal, second_length_axis, second_width_axis))


    # Generate ladder structures
    if VERBOSE:
        print(f"\n--------------------------------------------------")
        print(f"    Generating {n_ladder_structures} ladder-like structures:")
        print(f"--------------------------------------------------")
    for m in range(n_ladder_structures):

        print(f"    Ladder structure {m+1}")
        # Characterise the first main fault
        first_center = sample_new_center(centers_list, favoured_axis_1, favoured_axis_2, domain_scale=1, eps_center=eps_center)
        centers_list.append(first_center)
        first_length = 1.4 *np.random.uniform(L_min, L_max)
        first_width = np.random.uniform(W_min, W_max)
        first_normal = random_unit_vector()
        first_length_axis, first_width_axis = plane_main_axes_from_normal(first_normal)
        faults_info.append((first_center, first_length, first_width, first_normal, first_length_axis, first_width_axis))            
        if VERBOSE:
            print(f"            Main parallel fault")
            display_fault_properties((first_center, first_length, first_width, first_normal, first_length_axis, first_width_axis))

        # Characterise the second main fault
        dist_between_parallels = np.min(np.array([np.random.uniform(.8*first_width, 1*first_width), 0.6*first_length]))
        second_center = first_center + dist_between_parallels * np.array(first_normal)
        second_length = first_length
        second_width = first_width
        second_normal = first_normal
        second_length_axis = first_length_axis
        second_width_axis = first_width_axis
        faults_info.append((second_center, second_length, second_width, second_normal, second_length_axis, second_width_axis))

        if VERBOSE:
            print(f"            Second parallel fault")
            display_fault_properties((second_center, second_length, second_width, second_normal, second_length_axis, second_width_axis))

        # Characterise smaller, middle faults: 
        N_smaller_perpendicular_faults = 3
        dist_between_smaller_faults = first_length / (N_smaller_perpendicular_faults+2)
        mid_length = 1.3 *dist_between_parallels
        mid_width = dist_between_parallels / 4

        # They are at a slight angle to the normal: 
        angle = np.random.choice([-1, 1]) * np.random.uniform(10, 20, size=1)
        alpha = np.cos(np.radians(90 - angle))
        beta = np.sqrt(1-alpha**2)
        mid_length_axis = beta*np.array(first_normal)  + alpha*np.array(first_length_axis)
        mid_length_axis = mid_length_axis / np.linalg.norm(mid_length_axis)
        mid_width_axis = first_width_axis
        mid_normal = np.cross(mid_length_axis, mid_width_axis) 
    
        for n in range(N_smaller_perpendicular_faults):
            if n % 2 == 0:
                m = n // 2
            else:
                m = -((n + 1) // 2)
            mid_center = (first_center + second_center)/2 + m*dist_between_smaller_faults*first_length_axis

            faults_info.append((mid_center, mid_length, mid_width, mid_normal, mid_length_axis, mid_width_axis))
            if VERBOSE:
                print(f"\n              Smaller middle fault: fault number {n}")
                display_fault_properties((mid_center, mid_length, mid_width, mid_normal, mid_length_axis, mid_width_axis))
        
    return faults_info


def generate_points_from_plane(center, length, width, normal, length_axis, width_axis, 
                                distribution_mode='normal', domain_scale=1, density = 150):
    """
    Samples points with the specified density/m^2 along/around the plane given by input fault metadata.

    Returns:
    - fault_points: List of 3D points belonging to the given plane
    """
    #density = 150 # points per m^2
    area = length*width 
    num_points = int(np.round(density * area)) 

    # Sample points coordinates on the plane
    x = np.random.uniform(-length/2, length/2, num_points).astype(float)
    y = np.random.uniform(-width/2, width/2, num_points).astype(float)

    # Sample distances from the plane
    # Convention for faults: mean distance ~ width => 4*sigma=width
    if distribution_mode == 'normal':
        distribution_param = width/7
        dist = np.random.normal(loc=0, scale=distribution_param, size=num_points)

    # Compute fault points
    fault_points = center + x[:, None] * length_axis + y[:, None] * width_axis + dist[:, None] * normal

    return fault_points



def generate_point_list_from_metadata_list(fault_metadata_list, n_simple_faults, n_bent_faults, n_cross_faults, n_Y_faults, n_parallel_faults):
    """
    Generates sampled point cloud for every fault given its metadata. 

    Returns:
    - points_list: list of sampled points
    - labels: labels corresponding to point clusters (starting at 0, up to n_faults)
    """
    points_list = []
    labels = []
    k = 0
    label = 0
    for i, metadata in enumerate(fault_metadata_list):

        # Add simple faults
        if (i <= n_simple_faults -1) and (n_simple_faults != 0): 
            (center, length, width, normal, length_axis, width_axis) = metadata
            fault_points = generate_points_from_plane(center, length, width, normal, length_axis, width_axis)
            points_list.append(fault_points)

            # Add labels
            fault_labels = [label] * len(fault_points)
            labels.extend(fault_labels)


        # Add bent faults, make sure to merge main-bent points into one list
        elif (i>=n_simple_faults) and (i <= n_simple_faults + 2 * n_bent_faults -1) and (n_bent_faults != 0): 
            if k%2 == 0: # main points
                (main_center, main_length, main_width, main_normal, main_length_axis, main_width_axis) = metadata
                
                main_fault_points = generate_points_from_plane(main_center, main_length, main_width, main_normal, main_length_axis, main_width_axis)
                k+=1

            else:        # bent points
                (second_center, second_length, second_width, second_normal, second_length_axis, second_width_axis) = metadata
                second_fault_points = generate_points_from_plane(second_center, second_length, second_width, second_normal, second_length_axis, second_width_axis)

                # Concatenate main and bent points into one array, and add to list.
                points_list.append(np.concatenate((main_fault_points, second_fault_points)))
                k+=1

                # Add labels
                fault_labels = [label] * (len(main_fault_points) + len(second_fault_points))
                labels.extend(fault_labels)


        # Add cross faults (a cross consists of two faults, hence they are separate)
        elif   (i >= n_simple_faults + 2*n_bent_faults) and (i <= n_simple_faults + 2*n_bent_faults + 2*n_cross_faults -1) and (n_cross_faults != 0):
            (center, length, width, normal, length_axis, width_axis) = metadata
            fault_points = generate_points_from_plane(center, length, width, normal, length_axis, width_axis, density = 200)
            points_list.append(fault_points)

            # Add labels
            fault_labels = [label] * len(fault_points)
            labels.extend(fault_labels)


        # Add Y-shaped faults
        elif   (i >= n_simple_faults + 2*n_bent_faults + 2*n_cross_faults) and (i <= n_simple_faults + 2*n_bent_faults + 2*n_cross_faults + 2*n_Y_faults -1) and (n_Y_faults != 0):
            (center, length, width, normal, length_axis, width_axis) = metadata
            fault_points = generate_points_from_plane(center, length, width, normal, length_axis, width_axis)
            points_list.append(fault_points)

            # Add labels
            fault_labels = [label] * len(fault_points)
            labels.extend(fault_labels)
            

        # Add parallel faults
        elif    (i >= n_simple_faults + 2*n_bent_faults + 2*n_cross_faults + 2*n_Y_faults) and (i <= n_simple_faults + 2*n_bent_faults + 2*n_cross_faults + 2*n_Y_faults -1 + 2*n_parallel_faults -1) and (n_parallel_faults != 0):
            (center, length, width, normal, length_axis, width_axis) = metadata
            fault_points = generate_points_from_plane(center, length, width, normal, length_axis, width_axis)
            points_list.append(fault_points)

            # Add labels
            fault_labels = [label] * len(fault_points)
            labels.extend(fault_labels)
           

        # Add rest of faults that form the ladder structure
        else: 
            (center, length, width, normal, length_axis, width_axis) = metadata
            fault_points = generate_points_from_plane(center, length, width, normal, length_axis, width_axis)
            points_list.append(fault_points)

            # Add labels
            fault_labels = [label] * len(fault_points)
            labels.extend(fault_labels)
            
        label += 1
        
    return points_list, labels



def generate_dataset(n_simple_faults, n_bent_faults, n_cross_faults, n_Y_faults, n_parallel_faults, n_structures, VERBOSE=False):

    fault_metadata_list = generate_faults_metadata(n_simple_faults, n_bent_faults, n_cross_faults, n_Y_faults, n_parallel_faults, n_structures, eps_center=0.5, domain_scale=1, VERBOSE=False)

    fault_points_list, labels = generate_point_list_from_metadata_list(fault_metadata_list, n_simple_faults, n_bent_faults, n_cross_faults, n_Y_faults, n_parallel_faults)
    
    return fault_points_list, labels




#############       Visualization Functions        ##############


def multi_static_plot(fault_points_list, fault_metadata_list=None, show_plane=False, show_axes=False):
    """
    Visualize fault points in 3D and optionally plot the plane and axes.

    Args:
    - fault_points_list: List of fault points arrays (each array represents points for one fault)
    - fault_metadata_list: List of fault metadata tuples (each tuple contains metadata for one fault)
    - show_plane: Boolean indicating whether to plot the plane (default is False)
    - show_axes: Boolean indicating whether to plot the axes (default is False)
    """
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection='3d')

    # Plot event points for each fault
    for i, fault_points in enumerate(fault_points_list):
        ax.scatter(fault_points[:, 0], fault_points[:, 1], fault_points[:, 2], label=f'Fault {i+1}', s=5)

    # Plot the plane and axes if requested
    if show_plane or show_axes or fault_metadata_list != None:
        for metadata in fault_metadata_list:
            # Unpack plane data
            (center, length, width, normal, length_axis, width_axis) = metadata

            if show_plane:
                # Define vertices of the plane
                plane_vertices = center + np.array([-length/2, length/2, length/2, -length/2, -length/2])[:, None] * length_axis + np.array([-width/2, -width/2, width/2, width/2, -width/2])[:, None] * width_axis

                # Plot the plane
                ax.plot_trisurf(plane_vertices[:, 0], plane_vertices[:, 1], plane_vertices[:, 2], color='blue', alpha=0.2)

            if show_axes:
                # Plot length and width axes
                ax.plot([center[0], center[0] + length_axis[0]*length/2],
                        [center[1], center[1] + length_axis[1]*length/2],
                        [center[2], center[2] + length_axis[2]*length/2], color='red', linewidth=3)

                ax.plot([center[0], center[0] + width_axis[0]*width/2],
                        [center[1], center[1] + width_axis[1]*width/2],
                        [center[2], center[2] + width_axis[2]*width/2], color='green', linewidth=3)

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title('Synthetic Seismic Data - Fault Points')
    ax.legend()
    plt.show()


def multi_interactive_plot(fault_points_list, fault_metadata_list, show_plane=False, show_axes=False, domain_scale=1):
    """
    Visualize 3D points of multiple faults in an interactive plot using Plotly.
    """

    fig = go.Figure()

    x_range_largest = [-1, 1]
    y_range_largest = x_range_largest
    z_range_largest = x_range_largest
    
    for i, points in enumerate(fault_points_list):
        # Extract x, y, z coordinates from points
        x_coords = points[:, 0]
        y_coords = points[:, 1]
        z_coords = points[:, 2]

        # Plot the points
        fig.add_trace(go.Scatter3d(
            x=x_coords,
            y=y_coords,
            z=z_coords,
            mode='markers',
            marker=dict(
                size=3,
                opacity=0.8
            ),
            name=f'Fault {i+1}'
        ))

    if show_plane or show_axes or fault_metadata_list != None:
        j=0
        print(fault_metadata_list[0])
        for metadata in fault_metadata_list:
            # Unpack plane data:
            (center, length, width, normal, length_axis, width_axis) = metadata

            if show_plane:
                # Define vertices of the plane
                plane_vertices = center + np.array([-length/2, length/2, length/2, -length/2, -length/2])[:, None] * length_axis + np.array([-width/2, -width/2, width/2, width/2, -width/2])[:, None] * width_axis

                plane_trace = go.Mesh3d(x=plane_vertices[:, 0],
                                        y=plane_vertices[:, 1],
                                        z=plane_vertices[:, 2],
                                        i=[0, 1, 2, 0],
                                        j=[1, 2, 3, 4],
                                        k=[2, 3, 4, 0],
                                        opacity=0.2,
                                        color='blue',
                                        )
                fig.add_trace(plane_trace)

            if show_axes:
                # Plot length and width axes
                length_axis_trace = go.Scatter3d(x=[center[0], center[0] + length_axis[0]*length/2],
                                                y=[center[1], center[1] + length_axis[1]*length/2],
                                                z=[center[2], center[2] + length_axis[2]*length/2],
                                                mode='lines',
                                                line=dict(color='red', width=3),
                                                )

                width_axis_trace = go.Scatter3d(x=[center[0], center[0] + width_axis[0]*width/2],
                                                y=[center[1], center[1] + width_axis[1]*width/2],
                                                z=[center[2], center[2] + width_axis[2]*width/2],
                                                mode='lines',
                                                line=dict(color='green', width=3),
                                                )
                fig.add_trace(length_axis_trace)
                fig.add_trace(width_axis_trace)
            j+=1
            
    fig.update_layout(scene=dict(
        xaxis_title='X',
        yaxis_title='Y',
        zaxis_title='Z'),
        title='3D Visualization of Fault Points'
    )

    fig.show()


def plot_label_distribution(train_labels):
    train_labels = np.array(train_labels, dtype=int)
    train_counts = np.bincount(train_labels)
    total_train_samples = len(train_labels)
    num_labels = len(train_counts)

    # Plot the bar plot
    plt.figure(figsize=(8, 5))
    x = np.arange(num_labels)

    plt.bar(x, train_counts, edgecolor='black')

    # Add text labels on the bars
    for i, count in enumerate(train_counts):
        plt.text(x[i], count + 0.1, str(count), ha='center', fontsize=10)

    plt.ylabel('Number of Samples')
    plt.title("Label Distribution")
    plt.xticks(x, range(num_labels))
    plt.tight_layout()
    plt.show()
    # plt.savefig(os.path.join(output_dir, f"{savename}"))
    plt.close()



def main():

    # Example command: python synthetic_data_generation_script.py -v --visualise --n_simple_faults 2 --n_bent_faults 2 --n_cross_faults 1 --n_Y_faults 1 --n_parallel_faults 0 --n_structures 0

    # Create the argument parser
    parser = argparse.ArgumentParser(description='Description of your script.')

    # Define arguments
    parser.add_argument('-o', '--output', type=str, help='Path to the output file')
    parser.add_argument('-v', '--verbose', action='store_true', help='Enable verbose mode')
    parser.add_argument('--visualise', action='store_true', help='Visualise generated dataset')
    parser.add_argument('--n_simple_faults', type=int, default=1, help='Number of simple faults')
    parser.add_argument('--n_bent_faults', type=int, default=0, help='Number of bent faults')
    parser.add_argument('--n_cross_faults', type=int, default=0, help='Number of cross faults')
    parser.add_argument('--n_Y_faults', type=int, default=0, help='Number of Y faults')
    parser.add_argument('--n_parallel_faults', type=int, default=0, help='Number of parallel faults')
    parser.add_argument('--n_structures', type=int, default=0, help='Number of structures')

    # Parse the arguments
    args = parser.parse_args()

    # Access the arguments
    output_file = args.output
    verbose = args.verbose
    visualise = args.visualise
    n_simple_faults = args.n_simple_faults
    n_bent_faults = args.n_bent_faults
    n_cross_faults = args.n_cross_faults
    n_Y_faults = args.n_Y_faults
    n_parallel_faults = args.n_parallel_faults
    n_structures = args.n_structures

    # Implement the main logic of your script here
    if verbose:
        print("Verbose mode enabled")
        print(f"Output file: {output_file}")
        print(f"n_simple_faults: {n_simple_faults}")
        print(f"n_bent_faults: {n_bent_faults}")
        print(f"n_cross_faults: {n_cross_faults}")
        print(f"n_Y_faults: {n_Y_faults}")
        print(f"n_parallel_faults: {n_parallel_faults}")
        print(f"n_structures: {n_structures}")

    fault_points_list, labels = generate_dataset(
        n_simple_faults=n_simple_faults,
        n_bent_faults=n_bent_faults,
        n_cross_faults=n_cross_faults,
        n_Y_faults=n_Y_faults,
        n_parallel_faults=n_parallel_faults,
        n_structures=n_structures,
        VERBOSE=verbose
    )

    if visualise:
        multi_static_plot(fault_points_list, None, show_plane=False, show_axes=False)
        # Uncomment if you want to use interactive plotting
        #multi_interactive_plot(fault_points_list, None, show_plane=False, show_axes=False)

    if output_file:
        # Write fault_points_list and labels to the output CSV file
        with open(output_file, 'w', newline='') as csvfile:
            csvwriter = csv.writer(csvfile)

            # Write the header
            csvwriter.writerow(['x', 'y', 'z', 'labels'])
            pts = np.concatenate(fault_points_list, axis=0)
            # Write the data
            for fault_point, label in zip(pts, labels):
                csvwriter.writerow(fault_point.tolist()+[label])
    else: 
        return fault_points_list, labels



if __name__ == "__main__":
    np.random.seed(71) #for yfault and ladder
    np.random.seed(31) #for cross
    main()