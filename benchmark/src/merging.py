import numpy as np
import multiprocessing as mp
from functools import partial

from kernelparameters import KernelParameters
from kernel_fitting import get_minimum_bbox
from kernel_measures import *
from merging_utils import *

def is_candidate(
    r: int,
    c: int,
    r_bbox: np.ndarray,
    c_bbox: np.ndarray,
    r_bkg: bool,
    c_bkg: bool,
    X: np.ndarray,
    labels: np.ndarray,
    align_t: float):

    ''' 
    Checks whether a given pair of kernels is a potential candidate for merging

    Parameters
    -----------
    r,c : int
        The indices of the two kernels

    r_bbox, c_bbox: np.ndarray
        The bounding boxes of the two kernels

    r_bkg, c_bkg: bool
        Indicate whether the two kernels are background kernels

    X: np.ndarray
        The data points

    labels: np.ndarray
        The kernel assignemnt of each data point

    align_t: float
        The threshold to use for the alignment score. If the score of a pair is
        smaller, it will not be considered a candidate for merging.

    Returns
    --------
    delete: bool
        Indicates whether the pair should be deleted from the candidate list
    '''
    
    delete = r_bkg or c_bkg

    if not delete:
        delete = not have_overlap(r_bbox, c_bbox)

    if align_t > 0 and not delete:
        align_score = determine_alignment(X[labels == r], X[labels == c], 71)
        delete = align_score < align_t

    return delete


def get_merging_candidates(X: np.ndarray,
                           kernels: KernelParameters,
                           kernel_assign: np.ndarray,
                           candidates: np.ndarray = None,
                           align_t: float = 0):
    
    '''  
    Get the ids of all candidate kernel pairs to be considered for merging based on overlapping bounding boxes

    Parameters
    -----------
    X: np.ndarray
        The data points as an array of shape (n_samples, n_dim)
        
    kernels: KernelParameters
        The current kernel configuration

    kernel_assign: np.ndarray
        The kernel assignment of each data point in X
    
    candidates: np.ndarray | None
        If the candidate matrix is provided, use it to remove pairs with nans from the candidates.
        If it is None, determine the candidates from scratch, default = None

    align_t: float
        The threshold to use for the alignment score. If the score of a pair is
        smaller, it will not be considered a candidate for merging.
        Must be in [0,1], default = 0

    Returns
    -------
    rows, cols: np.ndarray
        The ids of the two kernels in each pair
    
    candidates: np.ndarray
        A matrix of shape (n_kernels, n_kernels) with float entries and np.nan entries indicating
        whether each kernel pair is a candidate or not
    '''

    if align_t > 1 or align_t < 0:
        raise ValueError(f'Alignment threshold must be in [0,1], was {align_t}')
    
    if candidates is None:
        candidates = np.zeros((kernels.get_n_kernels(),kernels.get_n_kernels()), dtype=np.float64)

    # get indices of all relevant kernel pairs
    rows,cols = np.mask_indices(len(candidates), lambda x,k: np.logical_and(np.triu(x,k),~np.isnan(candidates)),1)

    #print(f'Number of nans in gain {np.sum(np.isnan(candidates))/2}/{len(candidates)*(len(candidates)-1)/2}')
    _, __, ___, r_bbox, r_bkg = kernels.get_kernels(rows)
    _, __, ___, c_bbox, c_bkg = kernels.get_kernels(cols)

    # cutoff for parallel execution
    parallel_threshold = 1000

    if len(rows) > parallel_threshold:

        print('  Running in parallel')

        # use all but one available logical cpus
        n_cpus = mp.cpu_count()-1
        pool = mp.Pool(processes = n_cpus)

        del_idx = []
        pair_params = zip(rows, cols, r_bbox, c_bbox, r_bkg, c_bkg)

        # apply candidate function in parallel
        del_idx = np.array(pool.starmap(partial(is_candidate, X = X, labels = kernel_assign, align_t = align_t), pair_params, chunksize = parallel_threshold//2))

    else:
        # apply candidate function sequentially
        del_idx = np.array(list(map(partial(is_candidate, X = X, labels = kernel_assign, align_t = align_t), rows, cols, r_bbox, c_bbox, r_bkg, c_bkg)))
    
    # remove all non-candidates
    candidates[rows[del_idx],cols[del_idx]] = np.nan
    candidates[cols[del_idx],rows[del_idx]] = np.nan

    #print(f'Nans from bbox check: {sum(del_idx)}')
    #print(f'Number nans after bbox check: {np.sum(np.isnan(candidates))/2}')

    rows = rows[~del_idx]
    cols = cols[~del_idx]

    print(f'  {sum(~del_idx)}/{len(del_idx)} pairs have touching bbox')

    return rows, cols, candidates




def merge_single_pair(kernel_pair: KernelParameters,
                      keep_wgt: bool = True):
    '''
    Perform a Gaussian merge on a pair of kernels

    Parameters
    ----------
    kernel_pair: KernelParameters
        Contains the two kernels for which the merged kernel should be computed
    
    keep_wgt: bool
        Indicates whether to compute the weight of the merged kernel or set it to 1, default = True

    Returns
    --------
    mean: np.ndarray
        The mean vector of the merged kernel, has shape (1, n_dim)
    
    cov: np.ndarray
        The covariance matrix of the merged kernel, has shape (1, n_dim, n_dim)
    
    weight: np.ndarray
        The weight of the merged kernel, has shape (1,)
    
    bbox: np.ndarray
        The bounding box of the merged kernel, has shape(1, 2**n_dim, n_dim)
    
    is_bkg: np.ndarray
        Indicates whether the merged kernel is a background kernel, has shape (1,)

    '''
    if kernel_pair.get_n_kernels() != 2:
        raise ValueError('Kernel configuration is not a pair')

    # unpack the kernel parameters
    m1, c1, w1, b1, ib1 = kernel_pair.get_kernels(0)
    m2, c2, w2, b2, ib2 = kernel_pair.get_kernels(1)
    
    # calculate merged bbox and weight
    weight = w1+w2
    bbox, center = get_minimum_bbox(np.concatenate([b1,b2], axis=0))
    is_bkg = ib1 or ib2

    # calculate new mean and covariance
    if is_bkg:
        mean = center
        cov = np.cov(bbox, rowvar=False)*3.5
    
    else:
        mean = 1/weight * (w1 * m1 + w2 * m2)
        cov = (w1/weight)*(c1 + np.outer((m1-mean),(m1-mean))) + (w2/weight)*(c2 + np.outer((m2-mean),(m2-mean)))


    if not keep_wgt:
        weight = 1

    return np.array([mean]), np.array([cov]), np.array([weight]), np.array([bbox]), np.array([is_bkg])



def score_candidate(X: np.ndarray,
                    kernel_pair: KernelParameters,
                    mode: str = 'global',
                    bic_init: float = None,
                    tot_prob_init: np.ndarray = None,
                    n_kernels: int = None,
                    lab_cts: np.ndarray = None):

    ''' 
    Get the score and merged kernel for a specifc pair of kernels.

    Parameters
    ----------
    X: np.ndarray

    kernel_pair: KernelParameters
        Contains the two candidate kernels for which the merge score should be computed
    
    mode: str
        The mode to use for the score computation. Must be 'global' or 'local', default = 'global'

    bic_init: np.float | None
        The initial BIC score of all kernels without merging. Has to be provided when mode = 'global'
    
    tot_prob_init: np.ndarray | None
        The total probability under all kernels without merging, of shape (n_samples, 1). Has to be provided when mode = 'global'

    n_kernels: int | None
        The total number of kernels without merging. Has to be provided when mode = 'global'

    lab_cts: np.ndarray | None
        The number of points assigned to the individual kernels. Has to be provided when mode = 'local'


    Returns
    --------
    tot_prob_merged, tot_prob_separate: np.ndarray
        The cumulative probability of each point under the merged kernel and the kernel pair

    score: float
        The score for merging the candidate pair

    new_kern: KernelParameters
        The parameters of the merged kernel
    '''

    if mode == 'global' and (bic_init is None or n_kernels is None or tot_prob_init is None):
        raise ValueError('Must specify bic_init, tot_prob_init and n_kernels in global mode')
    
    if mode == 'local' and lab_cts is None:
        raise ValueError('Must specify lab_cts in local mode')

    
    # calculate parameters of the merged kernel
    new_kern = KernelParameters()
    new_kern.add_kernels(*merge_single_pair(kernel_pair, keep_wgt = mode == 'global'))

    if mode == 'local':
        # modify weight
        kernel_pair.modify_weight(np.array(lab_cts, dtype=np.float64)/len(X))


    # get probability and bic under the merged kernel
    tot_prob_merged = get_kernel_prob(X, new_kern)

    # get probability and bic under the two separate kernels
    tot_prob_separate = get_kernel_prob(X, kernel_pair)

    if mode == 'global':

        tot_prob_separate = np.sum(tot_prob_separate, axis=1, keepdims = True)
        tot_prob_merged = np.sum(tot_prob_merged, axis=1, keepdims = True)

        sum_tot_prob = tot_prob_init - tot_prob_separate + tot_prob_merged
        bic_merged = get_bic(sum_tot_prob, n_kernels = n_kernels-1)
        bic_separate = bic_init

    else:
        bic_separate = get_bic(tot_prob_separate)
        bic_merged = get_bic(tot_prob_merged)
            
    
    # calculate the information gain from merging
    score = bic_separate - bic_merged
    
    return tot_prob_merged, tot_prob_separate, score, new_kern



def merge_clusters(X: np.ndarray,
                   kernels: KernelParameters,
                   kernel_assign: np.ndarray,
                   init_prob: np.ndarray,
                   cutoff_align_check: int,
                   gain_mode: str = 'global',
                   align_t: float = 0):
    
    '''
    Merge clusters iteratively based on their BIC

    Parameters
    -----------
    X: np.ndarray
        The datapoints of shape (n_samples, n_features)

    kernels: KernelParameters
        The current Kernel configuration

    kernel_assign: np.ndarray
        The cluster assignment of each datapoint as and array of shape (n_samples,)

    init_prob: np.ndarray
        The initial probability of each datapoint under each kernel as an array
        of shape (n_samples, n_kernels)

    cutoff_align_check: int
        The number of clusters below which to consider the alignment for determining candidates

    gain_mode: str
        The mode to use for calculating the scores of candidate pairs, default = 'global'

    align_t: float
        The alignment threshold below which two clusters are not considered for merging.
        Must be in [0,1], default = 0

    
    Returns
    -------
    kernels: KernelParameters
        The kernel configuration after merging
    
    kernel_assign: np.ndarray
        The cluster assignment of each datapoint as and array of shape (n_samples,).
        Will be identical to input if gain_mode = 'global'

    '''

    modes = ['local', 'global']

    if gain_mode not in modes:
        raise ValueError(f'Unknown gain mode {gain_mode}. Must be one of {modes}.')

    bic_init = get_bic(init_prob)
    tot_prob_init = np.sum(init_prob, axis=1, keepdims=True)
    n_points = len(X)
    
    its = 0
    gain = None

    while True:

        print(f'  BIC: {bic_init}')

        # do not check alignment in the beginning
        if kernels.get_n_kernels() > cutoff_align_check:
            align_t_curr = 0
        else: 
            align_t_curr = align_t
        
        align_t_curr = align_t

        # determine kernel pairs with touching bboxes and proper alignment
        rows, cols, gain = get_merging_candidates(X, kernels, kernel_assign, gain, align_t_curr)
        
        if len(rows)==0:
            print('  No candidate pairs')
            break


        new_kernels = KernelParameters()
        p_merged = []
        p_separate = []
        merge_score = []


        for r, c in zip(rows,cols):

            # get the pair kernels
            old_kerns = KernelParameters()
            old_kerns.add_kernels(*kernels.get_kernels([r,c]))
                
            if gain_mode == 'local':
                
                # only consider contributions from points assigned to the kernel pair
                X_curr = X[np.logical_or(kernel_assign == r, kernel_assign == c)]

                cts = np.array([sum(kernel_assign == r), sum(kernel_assign == c)])
            
            else:
                X_curr = X
                cts = None

            # score the candidate pair
            prob_merg, prob_sep, score, new_kern = score_candidate(X_curr, old_kerns, gain_mode, bic_init, tot_prob_init, kernels.get_n_kernels(), cts)

            # save the candidate stats
            p_merged.append(prob_merg)
            p_separate.append(prob_sep)
            merge_score.append(score)

            # save the new kernel
            new_kernels.concatenate_parameters(new_kern)

        if np.nanmax(merge_score) <= 0:
            
            print('  No pairs with information gain')
            break

        merge_score = np.array(merge_score)
        neg_idx = merge_score <= 0

        print(f'  Removed {sum(neg_idx)}/{len(neg_idx)} candidate pairs with negative gain')
        
        # write the new scores
        gain[rows, cols] = merge_score
        gain[cols, rows] = merge_score

        gain[rows[neg_idx], cols[neg_idx]] = np.nan
        gain[cols[neg_idx], rows[neg_idx]] = np.nan

        print(f'Nans from negative check: {sum(neg_idx)}')
        print(f'Number nans after negative check: {np.sum(np.isnan(gain))/2}')
    

        # get unique pairs
        rows, cols, merge_score, idx_sort = get_disjoint_pairs(rows, cols, np.array(merge_score), gain)
        n_pairs = len(rows)


        if gain_mode == 'local':

            # merge clusters with global method
            p_merged = np.zeros((n_points, 1))
            p_separate = np.zeros((n_points, 1))
            
            for r,c in zip(rows, cols):

                # get the pair kernels
                old_kerns = KernelParameters()
                old_kerns.add_kernels(*kernels.get_kernels([r,c]))
                
                # calculate parameters of the merged kernel
                new_kern = KernelParameters()
                new_kern.add_kernels(*merge_single_pair(old_kerns, keep_wgt = True))

                # get cumulative probability under the merged kernel
                p_merged += np.sum(get_kernel_prob(X, new_kern), axis=1, keepdims=True)

                # get cumulative probability under the two separate kernels
                p_separate += np.sum(get_kernel_prob(X, old_kerns), axis=1, keepdims=True)
                
                new_kernels.concatenate_parameters(new_kern)
            
        else:

            # transform sort index to boolean mask to use on KernelParameters
            msk_keep = np.full(new_kernels.get_n_kernels(), False)
            msk_keep[idx_sort] = True

            # can reuse candidate stats from above
            new_kernels.delete_kernels(msk_keep)
            p_merged = sum(np.array(p_merged)[msk_keep])
            p_separate = sum(np.array(p_separate)[msk_keep])


        # update the initial probability and get the bic of all the merging events
        tot_prob_init = tot_prob_init - p_separate + p_merged
        bic_init = get_bic(tot_prob_init, n_kernels=kernels.get_n_kernels()-new_kernels.get_n_kernels())

        if np.isinf(bic_init):
            print('  Invalid BIC detected')
            break

        # modify kernel assignment
        if gain_mode == 'local':
            kernel_assign = merge_kernel_assignment(kernel_assign, rows, cols, kernels.get_n_kernels())

        # delete the now merged kernels
        idx_del = np.concatenate([rows, cols], axis=0)
        kernels.delete_kernels(idx_del)

        # add the new kernels
        kernels.concatenate_parameters(new_kernels)

        # adjust the gain matrix
        gain = np.delete(gain, idx_del, axis = 0)
        gain = np.delete(gain, idx_del, axis = 1)

        gain = np.concatenate([gain, np.zeros((len(gain), new_kernels.get_n_kernels()))], axis=1)
        gain = np.concatenate([gain, np.zeros((new_kernels.get_n_kernels(), gain.shape[1]))], axis=0)
        
        print(f'  Merged {n_pairs} pairs  >> {kernels.get_n_kernels()} kernels left')
        its+=1

    print(f'  Merging terminated after {its} iterations')    
    return kernels, kernel_assign
