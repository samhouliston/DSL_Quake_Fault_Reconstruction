import numpy as np
import warnings
from scipy.stats import multivariate_normal

from kernelparameters import KernelParameters
from inhull import inhull


def get_kernel_prob(X: np.ndarray,
                     kernels: KernelParameters
                     ):
    '''
    Determine the probability of each data point under each kernel

    Parameters
    ----------
    X: np.ndarray
        An array of observations. Must be of shape (n_samples, n_dimensions).

    kernels: KernelParameters
        The current kernel configuration
        

    Returns
    --------
    kernel_prob: np.array
        An array of shape (n_samples, n_kernels)
    '''

    kernel_prob = np.zeros((len(X), kernels.get_n_kernels()), dtype=np.float64)

    # calculate the the probability of each point under each kernel
    for i in range(kernels.get_n_kernels()):

        mean, cov, weight, bbox, is_bkg = kernels.get_kernels(i)

        if not is_bkg:
            
            kernel_prob[:,i] = multivariate_normal(mean=mean, cov=cov).pdf(X) * weight

        else:
            # determine points in background kernel
            eps = 1e-10*np.mean(np.abs(bbox))
            bkg_pts = inhull(X, bbox, eps)
            
            # calculate probability of being in the background
            kernel_prob[bkg_pts,i] = weight/np.prod(np.sqrt(np.linalg.eigvalsh(cov)), axis=0)

    return kernel_prob




def get_bic(kernel_prob: np.ndarray,
            n_kernels: int = None):
    '''
    Calculate the Bayesian Information Criterion of the dataset from the probability of each datapoint under each kernel

    Parameters
    ----------
    kernel_prob: np.array
        An array of shape (n_samples, n_kernels) that contains the probability of every sample under every kernel
    
    n_kernels: int
        Optional number of kernels to use. Can be specified if kernel_prob is already a cumulative probability, default = None
    
    Returns
    --------
    bic: float
        The BIC value of the current kernel configuration
    '''
    
    n_points = kernel_prob.shape[0]
    
    if n_kernels is None:
        n_kernels = kernel_prob.shape[1]

    # get cumulative probability per data point
    total_prob = np.sum(kernel_prob, axis = 1)

    res_msk = total_prob < np.finfo(np.float64).eps
    if sum(res_msk) > 0:
        total_prob[res_msk] = np.finfo(np.float64).eps
        warnings.warn('Reached probability density resolution limit')

    # calculate BIC
    bic = np.sum(-np.log(total_prob))+0.5*(10*n_kernels-1)*np.log(n_points)

    return bic