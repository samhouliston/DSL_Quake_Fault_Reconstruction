import numpy as np
import open3d as o3d

from kernelparameters import KernelParameters
from minboundbox import get_bbox
from kernel_measures import get_kernel_prob

USE_LIB = False  # whether to use library for bbox calculation


def get_minimum_bbox(X: np.ndarray)->tuple:

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

    if USE_LIB:
        # create a point cloud object from the data
        cloud = o3d.geometry.PointCloud()
        cloud.points = o3d.utility.Vector3dVector(X)

        # get the corners and center of the minimum bounding box
        bbox = cloud.get_minimal_oriented_bounding_box()

        corners = np.asarray(bbox.get_box_points())
        center = bbox.get_center()

    else:
        corners = get_bbox(X)
        center = np.mean(corners, axis = 0)

    return corners, center



def get_gaussian_bbox(mean: np.ndarray,
                      cov: np.ndarray,
                      n_var: float)->tuple:
    '''
    Calculate the minimum volume oriented bounding box for the points within sqrt(n_var) standard deviations of the mean of a 3D-Gaussian.

    Parameters
    ----------

    mean: np.array
        Mean vector of the Gaussian. Must be of shape (3,).
    
    cov: np.array
        Covariance matrix of the Gaussian. Must be of shape (3,3).
    
    n_var: float
        How many variances from the mean of the Gaussian to consider for the bounding box.

    Returns
    --------

    corners: np.array
        The 8 corner points of the bounding box as an array of shape (8,3)

    mean: np.array
        The center point of the bounding box as an array of shape (3,). Is equivalent to mean vector of the Gaussian.

    '''

    if mean.shape != (3,) or cov.shape != (3,3):
        raise ValueError('Gaussian must be 3-dimensional')

    evals, evecs = np.linalg.eigh(n_var*cov)

    # determine lengths of the bbox edges
    l_sides = 0.5*np.sqrt(evals) 
    
    # get an unrotated cuboid centered around 0 of the correct size
    corners = np.array([[-1, 1, 1, -1, -1, 1, 1, -1],
                        [1, 1, 1, 1, -1, -1, -1, -1],
                        [-1, -1, 1, 1, 1, 1, -1, -1]], dtype=np.float64)
    corners *= np.expand_dims(l_sides, axis=1)
    
    # rotate and shift the cuboid to the correct position and get the correct shape
    corners = (evecs @ corners).T + np.tile(mean, [8,1])

    return corners, mean 




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
    kernels: KernelParameters
        The kernel configuration after fitting Gaussian kernels to valid clusters

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
    bbox = np.full((n_clusters+int(fit_background), 8, 3), np.nan)
    is_bkg = np.full(n_clusters+int(fit_background), False)

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

        bbox[-1,:,:], mean[-1,:] = get_minimum_bbox(X)
        is_bkg[-1] = True

        # compute background covariance
        covar[-1,:] = np.cov(bbox[-1,:,:], rowvar=False)*3.5 

        # compute background weight
        weight[-1] = sum([l not in valid_clusters for l in cluster_labels])/n_points
        
    kernels = KernelParameters()
    kernels.add_kernels(mean, covar, weight, bbox, is_bkg)


    return kernels




def assign_to_kernel(X: np.array,
                     kernels: KernelParameters,
                     min_sz_cluster: int = 4,
                     refit_gauss: bool = False,
                     refit_bbox: bool = True
                     ):
    '''
    Perform a single step of expectation maximization to assign each data point to its kernel.
    Assigns each point to the most likely kernel and then updates the kernels.

    Parameters
    ----------
    X: np.ndarray
        An array of observations. Must be of shape (n_samples, n_dimensions).

    kernels: KernelParameters
        The current kernel configuration
    
    min_sz_cluster: int
        The threshold on the cluster size for a cluster to be considered valid, default = 4

    refit_gauss: bool
        Indicates whether the mean and covariance of non-background kernels should be refitted, default = False

    refit_bbox: bool
        Indicates whether the bounding box of non-background kernels should be refitted, default = True
        
    Returns
    --------
    kernels: KernelParameters
        The updated kernel configuration
    
    k_assign: np.ndarray
        An array of shape (n_samples,) which indicates the kernel assignment of every point in X

    kernel_prob: np.ndarray
        An array of shape (n_samples, n_kernels) that contains the probability of each datapoint under each kernel
    '''

    n_points = len(X)

    # assign each data point to most likely kernel
    kernel_prob = get_kernel_prob(X, kernels)
    k_assign = np.argmax(kernel_prob, axis=1)

    del_ker = []

    for i in range(kernels.get_n_kernels()):

        msk = k_assign == i
        cluster_sz = sum(msk)
        mean, cov, weight, bbox, is_bkg = kernels.get_kernels(i)

        if cluster_sz >= min_sz_cluster or (is_bkg and cluster_sz > 0):
            
            if not is_bkg:

                if refit_gauss:
                    # update kernel parameters
                    cov = np.cov(X[msk,:], rowvar=False)
                    mean = np.mean(X[msk,:], axis=0)

                if refit_bbox:
                    # bbox defined by the points in the cluster
                    bbox_pts, _ = get_minimum_bbox(X[msk,:])

                    # bbox defined by the Gaussian kernel
                    bbox_gauss, _ = get_gaussian_bbox(mean, cov, 12)

                    # bbox is minimum bbox of the union of all corner points
                    bbox, _ = get_minimum_bbox(np.concatenate((bbox_pts, bbox_gauss), axis=0))
            
            weight = cluster_sz / n_points
            kernels.modify_kernels(i, mean, cov, weight, bbox, is_bkg)

        else:
    
            del_ker.append(i)

    # delete non-background kernels with too few points
    if len(del_ker) > 0:
        print(f"  Removed {len(del_ker)}/{kernels.get_n_kernels()} kernels with < {min_sz_cluster} points")
        kernels.delete_kernels(del_ker)

    # update probabilities and kernel assignment
    kernel_prob = get_kernel_prob(X, kernels)
    k_assign = np.argmax(kernel_prob, axis=1)

    return kernels, k_assign, kernel_prob