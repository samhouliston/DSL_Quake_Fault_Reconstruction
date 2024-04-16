import numpy as np
from scipy.cluster.hierarchy import ward, fcluster


def cut_chunks(X: np.ndarray,
               n_chunks: int):
    '''
    Cut X into n_chunks chunks based on its agglomerative tree

    Parameters
    ----------
    X: np.ndarray
        An array of observations. Must be of shape (n_samples, n_dimensions).
    
    n_chunks: int
        The number of chunks that the data should be cut into


    Returns
    --------
    chunk_labs: np.ndarray
        An array of shape (n_samples,) that contains the chunk assignment of each datapoint
    '''

    if n_chunks == 1:
        chunk_labs = np.zeros(len(X))
    
    else:
        # compute the labels from the agglomerative link tree
        link_tree = ward(X)
        chunk_labs = fcluster(link_tree, link_tree[-n_chunks, 2], "distance")-1

    return chunk_labs



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
    prev = [0,0]

    for n_branch in range(1,len(link_tree)):

        # get cluster labels and determine the cluster sizes
        labels = fcluster(link_tree, link_tree[len(link_tree)-n_branch, 2], "distance")-1
        _ , counts = np.unique(labels, return_counts = True)
        
        n_clusters = sum(counts >= min_sz_cluster)

        # check whether capacity has improved
        if n_clusters >= capacity:
            capacity = n_clusters
            capacity_labels = labels
        
        # stop when the number of clusters drops
        elif n_clusters < prev[0] and n_clusters < prev[1]:
            break

        prev.pop()
        prev.insert(n_clusters, 0)

    return capacity_labels

