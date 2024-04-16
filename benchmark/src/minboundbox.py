import numpy as np
import itertools
from scipy.spatial import ConvexHull

PREC = 10   # precision for rounding values



def euler_angles(vec1, vec2, normal):

    # calculate beta
    beta = np.arcsin(np.sign(normal[:,0])*np.minimum(np.ones(len(normal)), np.abs(normal[:,0])))

    alpha = np.zeros_like(beta)
    gamma = np.zeros_like(beta)

    msk = normal[:,0] == 1

    # calculate alpha
    
    if sum(msk) > 0:
        alpha[msk] = np.arcsin(np.sign(vec2[msk,2])* \
            np.minimum(np.ones(sum(msk)), np.abs( vec2[msk,2]) ))
        
    if sum(~msk) > 0:
        alpha[~msk] = np.arccos( np.sign( normal[~msk, 2]/ np.cos(beta[~msk]))* \
            np.minimum(np.ones(sum(~msk)), np.abs( normal[~msk,2] / np.cos(beta[~msk]))))
        
        alpha_msk = np.sign(normal[:,1]) != np.sign(-np.sin(alpha) * np.cos(beta))
        alpha[np.logical_and(~msk, alpha_msk)] *= -1

        # calculate gamma
        singamma = -vec2[~msk,0] / np.cos(beta[~msk])
        gamma_msk = vec1[~msk,0] >= 0

        gamma[gamma_msk] = np.arcsin(np.sign(singamma[gamma_msk]) * \
                            np.minimum(np.ones(sum(gamma_msk)), np.abs(singamma[gamma_msk])))
        gamma[~gamma_msk] = - np.pi - np.arcsin(np.sign(singamma[~gamma_msk]) * \
                            np.minimum(np.ones(sum(~gamma_msk)), np.abs(singamma[~gamma_msk])))
        

    return alpha, beta, gamma




def minrect(pts):

    # get the hull vertices
    hull = ConvexHull(pts)
    vertices = np.append(hull.vertices, hull.vertices[0])

    pts = pts[vertices]

    rotation = lambda theta: np.array([
        [np.cos(theta), np.sin(theta)],
        [-np.sin(theta), np.cos(theta)]
    ])

    x = pts[:,0]
    y = pts[:,1]

    x_diff = x[1:]-x[:-1]
    y_diff = y[1:]-y[:-1]

    x_diff[np.abs(x_diff) < np.finfo(np.float64).eps * 100] = np.finfo(np.float64).eps * 100
    y_diff[np.abs(y_diff) < np.finfo(np.float64).eps * 100] = np.finfo(np.float64).eps * 100


    angle = np.arctan2(y_diff, x_diff)
    angle = np.unique(np.mod(angle,np.pi/2))


    global_min = np.inf
    global_rot = None

    # find the optimal rotation
    for a in angle:
        
        rot = rotation(-a)
        pts_rot = pts @ rot

        min_val = np.min(pts_rot, axis=0)
        max_val = np.max(pts_rot, axis=0)

        A = np.prod(max_val-min_val)
        A = np.round(A, PREC)

        if A < global_min:
            global_min = A
            global_rot = rot

    
    return global_rot




def checkbox(alpha, beta, gamma, X, min_vol, minmax, min_rot):

    # construct the rotation matrix from the euler angles
    rotation = np.array([
        [np.cos(beta)*np.cos(gamma), -np.cos(beta)*np.sin(gamma), np.sin(beta)],
        [np.sin(alpha)*np.sin(beta)*np.cos(gamma)+np.cos(alpha)*np.sin(gamma), 
         -np.sin(alpha)*np.sin(beta)*np.sin(gamma)+np.cos(alpha)*np.cos(gamma), -np.sin(alpha)*np.cos(beta)],
        [-np.cos(alpha)*np.sin(beta)*np.cos(gamma)+np.sin(alpha)*np.sin(gamma), 
         np.cos(alpha)*np.sin(beta)*np.sin(gamma)+np.sin(alpha)*np.cos(gamma), np.cos(alpha)*np.cos(beta)]
    ])

    # rotate the face into the xy-plane
    X_r = X @ rotation
    rot_opt = minrect(X_r[:,:2])

    rot_opt = np.concatenate([rot_opt, np.zeros((2,1))], axis = 1)
    rot_opt = np.concatenate([rot_opt, [[0,0,1]]], axis = 0)

    rotation = rotation @ rot_opt

    X_r = X @ rotation

    curr_min = np.min(X_r, axis = 0)
    curr_max = np.max(X_r, axis = 0)

    vol = np.round(np.prod(curr_max - curr_min), PREC)

    if vol < min_vol:

        min_vol = vol
        minmax = np.array([curr_min, curr_max])
        min_rot = rotation

    return min_vol, minmax, min_rot



def get_bbox(X: np.ndarray):

    X = np.round(X, PREC)

    hull = ConvexHull(X)

    simplices = np.sort(hull.simplices, axis = 1)

    # determine the 2 edges of every hull simplex
    edge1 = X[simplices[:,1]] - X[simplices[:,0]]
    edge1 /= np.repeat(np.sqrt(np.sum(edge1**2, axis=1, keepdims=True)), X.shape[1], axis=1)


    edge2 = X[simplices[:,2]] - X[simplices[:,1]]
    edge2 -= np.repeat(np.sum(edge1*edge2, keepdims=True, axis=1), X.shape[1], axis=1) * edge1
    edge2 /= np.repeat(np.sqrt(np.sum(edge2**2, axis=1, keepdims=True)), X.shape[1], axis=1)
    

    # get the normal vectors of the simplices
    normals = np.cross(edge1, edge2)


    alpha, beta, gamma = euler_angles(edge1, edge2, normals)

    min_vol = np.inf
    minmax = None
    min_rot = None

    # check all faces for the one with the smallest volume box
    for i in range(hull.nsimplex):

        min_vol, minmax, min_rot = checkbox(alpha[i], beta[i], gamma[i], X[hull.vertices], min_vol, minmax, min_rot)
    
    corner_pts = np.array([[minmax[x,0], minmax[y,1], minmax[z,2]] for x,y,z in itertools.product([0, 1], repeat=3)])
    corner_pts = corner_pts @ np.linalg.inv(min_rot)

    return corner_pts