import numpy as np 
import scipy as sc 
from scipy.cluster.hierarchy import ward, fcluster
import open3d as o3d


def get_capacity_clusters(X: np.array, 
                          min_sz_cluster: int = 4, 
                          min_n_merges: int = 4
                          )->np.array:
    '''
    Get the cluster assignment with the largest number of valid clusters based on ward linkage.
    Iteratively cuts the tree from the bottom until the number of valid clusters does not increase anymore.
    
    Parameters
    -----------
    X: np.array
        The data to compute the cluster assignment for as an array of size (n_samples x n_dimensions)
    
    min_sz_cluster: int
        The threshold on the cluster size for a cluster to be considered valid, default = 4

    min_n_merges: int
        The number of cluster merging steps that can be skipped by the algorithm, default = 4


    Returns
    --------
    capacity_labels: np.array
        The cluster labels as an array of length n_samples
    '''

    min_n_merges = max(min_sz_cluster, min_n_merges)

    link_tree = ward(X)

    capacity = 0
    capacity_labels = np.zeros(len(X))

    for n_merges in range(min_n_merges,len(X)):

        # get cluster labels and determine the cluster sizes
        labels = fcluster(link_tree, link_tree[n_merges-1, 2], "distance")
        uq_labs , counts = np.unique(labels, return_counts = True)
        
        n_clusters = sum(counts >= min_sz_cluster)

        # check whether capacity has improved
        if n_clusters >= capacity:
            capacity = n_clusters
            capacity_labels = labels
        
        # stop when all points are included in valid clusters
        elif n_clusters == len(uq_labs):
            break


    return capacity_labels


def get_minimum_bbox(X: np.array)->tuple:

    '''
    Calculate the minimum volume oriented bounding box for the points in X

    Parameters
    -----------
    X: np.array
        An array containing a point cloud of observations. Has to be of shape (n_samples, 3)

    Returns
    --------
    corners: np.array
        The 8 corner points of the bounding box as an array of shape (8,3)

    center: np.array
        The center point of the bounding box as an array of shape (3,)
    
    '''

    if X.shape[1] != 3:
        raise ValueError(f'Data has to be 3-dimensional, was {X.shape[1]}-dimensional')

    # create a point cloud object from the data
    cloud = o3d.geometry.PointCloud()
    cloud.points = o3d.utility.Vector3dVector(X)

    # get the corners and center of the minimum bounding box
    bbox = cloud.get_minimal_oriented_bounding_box()

    corners = np.asarray(bbox.get_box_points())
    center = bbox.get_center()

    return corners, center




def fit_gaussian_kernels(X: np.array,
                         cluster_labels: np.array,
                         min_sz_cluster: int = 4
                         )->tuple:
    '''
    Fit a Gaussian kernel to every valid cluster. Fits mean, covariance and weight for every kernel.
    If there are points that are not in any valid cluster, fit a uniform background kernel.

    Parameters
    -----------
    X: np.array
        An array of observations. Must be of shape (n_samples, n_dimensions).
    
    cluster_labels: np.array
        The cluster assignment of every observation. Must be of length n_samples
    
    min_sz_cluster: int
        The threshold on the cluster size for a cluster to be considered valid, default = 4

    Returns
    -------
    mean: np.array
        Array of shape (n_kernels, n_dimensions) that contains the covariance matrices for every kernel. 
        The last entry corresponds to the background kernel if there exist background points.

    covar: np.array
        Array of shape (n_kernels, n_dimensions, n_dimensions) that contains the covariance matrices for every kernel. 
        The last entry corresponds to the background kernel if there exist background points.
    
    weight: np.array
        Array of shape (n_kernels,) that contains the weight of every kernel.
        The last entry corresponds to the background kernel if there exist background points.

    bbox: None | np.array
        An array of shape (8,3) that contains the corners of the bounding box for the background kernel.
        None if there are no background points.

    '''

    if len(X) != len(cluster_labels):
        raise ValueError(f'Number of datapoints {len(X)} does not match number of labels {len(cluster_labels)}')
    
    # determine dataset parameters
    n_dim = X.shape[1]
    n_points = X.shape[0]
    uq_clusters, cluster_szs = np.unique(cluster_labels, return_counts = True) # all clusters and sizes
    n_clusters = sum(cluster_szs >= min_sz_cluster) # number of valid clusters
    valid_clusters = uq_clusters[cluster_szs >= min_sz_cluster] # labels of valid clusters
    fit_background = n_clusters < len(uq_clusters)

    mean = np.zeros((n_clusters+int(fit_background), n_dim))
    covar = np.repeat([np.eye(n_dim)], n_clusters+int(fit_background), axis=0)
    weight = np.zeros(n_clusters+int(fit_background))
    bbox = None

    for i in range(n_clusters):

        id = valid_clusters[i]
        X_curr = X[cluster_labels == id] 

        # compute cluster mean
        mean[i,:] = np.mean(X_curr, axis=0)

        # compute cluster covariance
        covar[i,:,:] = np.cov(X_curr, rowvar=False)

        # compute cluster weight
        weight[i] = len(X_curr)/n_points
    
    if fit_background and n_dim == 3:
        # fit the background kernel

        X_bkg = X[[l not in valid_clusters for l in cluster_labels]]
        bbox, center = get_minimum_bbox(X_bkg)

        # set background mean
        mean[-1,:] = center

        # compute background covariance
        #FIXME: sqrt(12)*stddev in paper but in the implementation it's sqrt(12)*variance?
        covar[-1,:] = np.cov(bbox, rowvar=False)*3.5 

        # compute background weight
        weight[-1] = len(X_bkg)/n_points
    

    return mean, covar, weight, bbox

        


