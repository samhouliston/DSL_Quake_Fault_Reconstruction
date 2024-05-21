import numpy as np
from scipy.io import loadmat
from sklearn.metrics import adjusted_rand_score

from kernelparameters import KernelParameters
from kernel_fitting import *
from capacity import *
from merging import *

import cProfile


def run_fault_reconstruction(X: np.ndarray,
                       min_sz_cluster: int,
                       n_chunks: int = 1,
                       gain_mode: str = 'global',
                       align_t: float = 0,
                       pct_align_check: float = 1,
                       margin_scale: float = 0,
                       min_n_clusters: int = 0
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

    align_t: float
        The alignment threshold below which two clusters are not merged.
        Must be in [0,1], default = 0

    pct_align_check: float
        Below which fraction of the capacity cluster number to start checking for the alignment
    
    min_n_clusters: int
        The minimum number of clusters to merge into, excluding the background, default = 0

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

        capacity = sum(~kernels.get("ib"))
        print(f'  Fit {capacity} Gaussian and {sum(kernels.get("ib"))} background kernels')

        # assign points and update kernels with EM
        kernels, cluster_labs, kernel_prob = assign_to_kernel(X[msk], kernels, min_sz_cluster, refit_gauss = False)

        # run the kernel merging algorithm
        kernels, cluster_labs = merge_clusters(X[msk], kernels, cluster_labs, kernel_prob, int(capacity*pct_align_check), gain_mode, align_t, margin_scale, min_n_clusters)

        # reassign the points with EM
        kernels, cluster_labs, kernel_prob = assign_to_kernel(X[msk], kernels, min_sz_cluster, refit_gauss = False)

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
    all_kernels, all_labels = merge_clusters(X, all_kernels, all_labels, kernel_prob, len(X), gain_mode, align_t, 0, min_n_clusters)
    
    # reassign the points with EM
    all_kernels, all_labels, kernel_prob = assign_to_kernel(X, all_kernels, min_sz_cluster, refit_gauss = False)

    print(f'Final number of kernels: {all_kernels.get_n_kernels()}')

    return all_kernels, all_labels



def main():

    run_case = 'Landers'
    
    if run_case == 'Landers':
        path = '../data/landers.mat'
        X = loadmat(path)['xyz_mat']

    if run_case == 'Synthetic':
        path = '../data/synthetics_31.npy'
        ground_truth = np.load(path)

        X = np.concatenate(ground_truth, axis = 0)
        
    kernels, labels = run_fault_reconstruction(X, min_sz_cluster = 4, n_chunks = 1, align_t = 0)


if __name__ == '__main__':
    main()