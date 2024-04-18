import numpy as np
from sklearn.decomposition import PCA

from inhull import inhull


def have_overlap(bbox1: np.ndarray,
                 bbox2: np.ndarray):
    '''
    Check whether two bounding boxes overlap.

    Parameters
    ----------
    bbox1: np.ndarray
        The corner points of the first bounding box as an array of shape (n_points, n_dimensions)
    
    bbox1: np.ndarray
        The corner points of the second bounding box as an array of shape (n_points, n_dimensions)
    
    Returns
    --------
    overlap: bool
        Indicates whether the two bounding boxes overlap

    '''
    
    overlap = np.any(inhull(bbox1, bbox2, 1.e-13*np.mean(np.abs(bbox2)))) or np.any(inhull(bbox2, bbox1, 1.e-13*np.mean(np.abs(bbox1))))

    return overlap




def get_disjoint_pairs(rows, cols, score, gain):

    '''
    Get the disjoint pairs from a set of row-column index pairs.
    The first occurrence of an index is decided based on descending score.
    All pairs that have positive gain with the result pairs are removed.

    Parameters
    -----------
    rows: np.ndarray
        An array of row indices
    
    cols: np.ndarray
        An array of column indices in the same order as rows
    
    score: np.ndarray
        An array of scores that determine the precedence of a pair in the same order as rows.
        Higher score is better and pairs with non-positive scores are cut.

    gain: np.ndarray
        A matrix of size (n_kernels, n_kernels) that contains all available merge scores of all kernels
    
    Returns
    -------
    rows, columns, score: np.array
        The row and column indices and scores of the unique pairs sorted by descending score
    
    idx_sort: np.array
        An index that can be used to extract the matching elements of an array with the same order
        as the input rows or columns
    '''

    # sort the pairs by descending score and remove negatives    
    idx_sort = np.argsort(score)[::-1]
    score = score[idx_sort]
    idx_imp = score > 0
    score = score[idx_imp]
    rows = rows[idx_sort][idx_imp]
    cols = cols[idx_sort][idx_imp]
    idx_sort = idx_sort[idx_imp]

    unique_rows = []
    unique_cols = []
    unique_score = []
    unique_idx = []

    while len(rows) > 0:
        
        # get the next pair with highest score
        r = rows[0]
        c = cols[0]

        unique_rows.append(r)
        unique_cols.append(c)
        unique_score.append(score[0])
        unique_idx.append(idx_sort[0])
        
        # find all kernels that could potentially be merged with the current pair
        curr_gain = gain[[r,c], :]
        curr_gain[np.isnan(curr_gain)] = -1
        exclude = np.argmax(curr_gain, axis = 0)
        exclude = np.array([idx for idx, i in enumerate(exclude) if gain[(r,c)[i], idx] > 0])
        
        for i in exclude:
            
            # remove all pairs that are in exclude
            if i in rows or i in cols:

                to_remove = np.logical_and(rows != i, cols != i)

                rows = rows[to_remove]
                cols = cols[to_remove]
                score = score[to_remove]
                idx_sort = idx_sort[to_remove]


    return np.array(unique_rows), np.array(unique_cols), np.array(unique_score), np.array(unique_idx)




def merge_kernel_assignment(kernel_assign: np.ndarray,
                            kernel1: np.ndarray, 
                            kernel2: np.ndarray,
                            n_old_kernels: int = None):
    
    '''
    Map the old kernel labels to the new kernel labels.
    Kernel labels coincide with position of the kernel in the KernelParameters object.

    Parameters
    -----------
    kernel_assign: np.ndarray
        The old kernel assignment of every data point. Has shape (n_samples,)
    
    kernel1, kernel2: np.ndarray
        Contain the labels of the first and second kernel of every merged pair in the same order

    n_old_kernels: int | None
        Manually specify the number of kernels before the merging. Needed if there are some kernels without any points assigned to them.
        If None, the number of kernels is inferred from the kernel assignment.


    Returns
    --------
    kernel_assign: np.ndarray
        The new kernel assignment of ever data point

    '''
    if n_old_kernels is None:
        n_old_kernels = len(np.unique(kernel_assign))

    # create a mapping from old kernel labels to new kernel labels
    old2new = np.arange(n_old_kernels)
        
    msk_keep = np.full(len(old2new), True)
    msk_keep[np.concatenate([kernel1, kernel2])] = False

    # add new labels of the old kernels 
    old2new = old2new[msk_keep]
    old2new = dict(zip(old2new, range(len(old2new))))

    n_old = len(old2new)

    # add new labels of the merged kernels
    for i in range(len(kernel1)):

        old2new[kernel1[i]] = n_old+i
        old2new[kernel2[i]] = n_old+i

    kernel_assign = np.array([old2new[l] for l in kernel_assign])

    return kernel_assign


def determine_alignment(X1: np.ndarray,
                        X2: np.ndarray,
                        seed: int
                        ):
    
    '''
    Determine the alignment of two point clouds based on their two leading principal components

    Parameters
    ----------
    X1, X2: np.ndarray
        The two point clouds as arrays of shape (n_samples, n_dim)
    
    seed: int
        The random seed

    Returns
    --------
    align_score: float
        The alignment score of the point clouds. Is calculated as the sum of 
        the dot products of their normalized leading principal components

    '''

    if X1.shape[1] != X2.shape[1]:
        raise ValueError('Point clouds must have the same number of features')

    # Determine the leading 2 principal components of each cloud
    pca1 = PCA(n_components = 2, 
              random_state = seed)
    pca2 = PCA(n_components = 2, 
              random_state = seed)
    
    comps1 = pca1.fit(X1).components_
    comps2 = pca2.fit(X2).components_

    # normalize the components
    comps1[0] /= np.linalg.norm(comps1[0])
    comps1[1] /= np.linalg.norm(comps1[1])
    comps2[0] /= np.linalg.norm(comps2[0])
    comps2[1] /= np.linalg.norm(comps2[1])

    # calculate alignment of the components
    align_score = sum(np.abs(np.sum(comps1 * comps2, axis = 1)))/2

    return align_score
