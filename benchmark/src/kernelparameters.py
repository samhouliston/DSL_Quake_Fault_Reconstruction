import numpy as np

class KernelParameters:
    '''
    A class to store the kernel parameters of all kernels in the fault network

    mean: np.array
        Array of shape (n_kernels, n_dimensions) that contains the mean for every kernel. 

    cov: np.array
        Array of shape (n_kernels, n_dimensions, n_dimensions) that contains the covariance matrices for every kernel. 
    
    weight: np.array
        Array of shape (n_kernels,) that contains the weight of every kernel.

    bbox: np.array
        An array of shape (n_kernels,8,3) that contains the corners of the bounding box.
    
    is_bkg: np.array
        A boolean array of shape(n_kernels,) that indicates which kernel is a background kernel
    '''

    def __init__(self, n_dim: int = 3):

        self.n_kernels = 0
        self.n_dim = n_dim
        
        self.mean = np.zeros((0,n_dim), dtype=np.float64)
        self.cov = np.zeros((0,n_dim,n_dim), dtype=np.float64)
        self.weight = np.zeros((0,), dtype=np.float64)
        self.bbox = np.zeros((0,2**n_dim,n_dim), dtype=np.float64)
        self.is_bkg = np.zeros((0,), dtype=bool)
    
    def get_dim(self):
        return self.n_dim
    
    def get_n_kernels(self):
        return self.n_kernels
    
    def add_kernels(self, mean, cov, weight, bbox, is_bkg):
        '''Add new kernels to the kernel configuration'''

        self.mean = np.concatenate([self.mean, mean], axis=0)
        self.cov = np.concatenate([self.cov, cov], axis=0)
        self.weight = np.concatenate([self.weight, weight], axis=0)
        self.bbox = np.concatenate([self.bbox, bbox], axis=0)
        self.is_bkg = np.concatenate([self.is_bkg, is_bkg], axis=0)

        self.n_kernels += len(weight)

    def concatenate_parameters(self, kernelp):
        '''Add new kernels from another KernelParameters object'''

        self.add_kernels(kernelp.mean, kernelp.cov, kernelp.weight, kernelp.bbox, kernelp.is_bkg)
    
    def modify_kernels(self, kernel_idx, mean, cov, weight, bbox, is_bkg):
        '''
        Modify a single or multiple existing kernels in the kernel configuration.
        Make sure that the type of kernel_idx and the shape of the other arguments match.
        '''

        self.mean[kernel_idx] = mean
        self.cov[kernel_idx] = cov
        self.weight[kernel_idx] = weight
        self.bbox[kernel_idx] = bbox
        self.is_bkg[kernel_idx] = is_bkg
    
    def modify_weight(self, weight):
        '''
        Modify the weights of all kernels in the kernel configuration
        '''

        if self.weight.shape != weight.shape:
            raise ValueError('Input weight matrix must have same shape as current weight matrix')

        self.weight = weight

    def get(self, field: str):
        '''
        Get a specific field from the kernel parameters.
        '''

        if field == 'm':
            return self.mean
        elif field == 'c':
            return self.cov
        elif field == 'w':
            return self.weight
        elif field == 'b':
            return self.bbox
        elif field == 'ib':
            return self.is_bkg
        else:
            raise ValueError(f'Unknown field {field}')
    
    def delete_kernels(self, kernel_idx):
        '''
        Remove a single or multiple kernels from the kernel configuration.

        Parameters
        -----------
        kernel_idx: int | array_like
            Integer, an integer index or a boolean mask of length n_kernels that indicates which kernels to keep

        '''

        if hasattr(kernel_idx, "__len__") and isinstance(kernel_idx[0], (bool, np.bool_)):
            self.mean = self.mean[kernel_idx]
            self.cov = self.cov[kernel_idx]
            self.weight = self.weight[kernel_idx]
            self.bbox = self.bbox[kernel_idx]
            self.is_bkg = self.is_bkg[kernel_idx]

            self.n_kernels = sum(kernel_idx)

        else:
            self.mean = np.delete(self.mean, kernel_idx, axis=0)
            self.cov = np.delete(self.cov, kernel_idx, axis=0)
            self.weight = np.delete(self.weight, kernel_idx, axis=0)
            self.bbox = np.delete(self.bbox, kernel_idx, axis=0)
            self.is_bkg = np.delete(self.is_bkg, kernel_idx, axis=0)

            self.n_kernels -= len(kernel_idx) if hasattr(kernel_idx, "__len__") else 1
    
    def get_kernels(self, kernel_idx = None):
        '''
        Extract the kernel parameters of a single or multiple kernels.
        If kernel_idx = None, return all kernels. 
        If kernel_idx is an integer, the outer dimension of the return values is flattened.
        '''

        if kernel_idx is None:
            kernel_idx = range(self.n_kernels)

        return self.mean[kernel_idx], self.cov[kernel_idx], self.weight[kernel_idx], self.bbox[kernel_idx], self.is_bkg[kernel_idx]