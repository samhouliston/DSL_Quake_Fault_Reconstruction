import numpy as np
from scipy.io import loadmat

from kernelparameters import KernelParameters
from kernel_fitting import *
from capacity import *
from merging import *


def run_fault_reconstruction(X: np.ndarray,
                       min_sz_cluster: int,
                       n_chunks: int = 1,
                       gain_mode: str = 'global'
                       ):
    '''
    Run the fault reconstruction algorithm (Kamer 2020)
    
    Parameters
    -----------
    X: np.ndarray
        The data to compute the cluster assignment for as an array of size (n_samples x n_dimensions)
    
    min_sz_cluster: int
        The threshold on the cluster size for a cluster to be considered valid, default = 4

    n_chunks: int
        Number of chunks that the data is cut into before performing the algorithm, default = 1
    
    gain_mode: str
        The method to use for determining the information gain of a merging event, default = 'global'

    Returns
    --------
    all_kernels: np.ndarray
        The cluster assignment of every data point
    
    all_labels: np.ndarray
        The cluster assignment of every data point
    '''

    n_points = len(X)

    # get the chunk partition
    chunk_labs = cut_chunks(X, n_chunks)

    all_kernels = KernelParameters()
    all_labels = np.empty(0, dtype=int)
    lab_offset = 0

    for chunk_id in range(n_chunks):

        msk = chunk_labs == chunk_id

        print(f'Processing chunk {chunk_id+1}/{n_chunks} with {sum(msk)} points')

        # get the capacity clusters based on the agglomerative tree
        capacity_labs = get_capacity_clusters(X[msk], min_sz_cluster)

        # fit the kernels
        kernels = fit_gaussian_kernels(X[msk], capacity_labs, min_sz_cluster)

        print(f'  Fit {sum(~kernels.get("ib"))} Gaussian and {sum(kernels.get("ib"))} background kernels')

        # assign points and update kernels with EM
        kernels, cluster_labs, kernel_prob = assign_to_kernel(X[msk], kernels, min_sz_cluster, refit_gauss = False)

        # run the kernel merging algorithm
        kernels, cluster_labs = merge_clusters(X[msk], kernels, cluster_labs, kernel_prob, gain_mode)

        # reassign the points with EM
        kernels, cluster_labs, kernel_prob = assign_to_kernel(X[msk], kernels, min_sz_cluster)

        print(f'  Merged capacity kernels into {kernels.get_n_kernels()} kernels')

        # shift and persist the cluster labels of the current chunk
        cluster_labs = cluster_labs + lab_offset
        lab_offset += kernels.get_n_kernels()
        all_labels = np.concatenate([all_labels, cluster_labs], axis=0)

        # scale the kernel weights
        kernels.modify_weight(kernels.get('w')*sum(msk)/n_points)

        # persist the kernels of the current chunk
        all_kernels.concatenate_parameters(kernels)
    

    # merge the kernels of all chunks
    print('Combine results of all chunks')
    kernel_prob = get_kernel_prob(X, all_kernels)
    all_kernels, all_labels = merge_clusters(X, all_kernels, all_labels, kernel_prob, gain_mode)
    print(f'Final number of kernels: {all_kernels.get_n_kernels()}')

    return all_kernels, all_labels



def main():

    
    path = '../data/landers.mat'
    X = loadmat(path)['xyz_mat']
    kernels, labels = run_fault_reconstruction(X, min_sz_cluster = 4)

if __name__ == '__main__':
    main()