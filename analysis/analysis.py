import numpy as np
import random as ran
import matplotlib.pyplot as plt
from tqdm import tqdm
from numba import njit

# Numba does not support arguments in np.var(), np.std()
# Therefore, var() and std() are equivalent to np.var(ddof=1) and np.std(ddof=1)

@njit
def var(data):
    n = data.size
    return data.var() * n/(n-1)

@njit
def std(data):
    return np.sqrt(var(data))

@njit
def halve(data):
    '''halve(numpy.ndarray 1D x): (x[i]) -> (y[i] = (x[2i] + x[2i+1])/2)
    y.size = x.size // 2
    If x.size is odd, excess elements are discarded at the end of the array.'''
    s = data.size
    w = np.copy(data)
    if (s % 2) != 0:
        w = w[1:]
        s -= 1
    s = s // 2
    return w.reshape(s, 2).sum(axis=1) / 2

@njit
def blocking(data):
    '''Blocking analysis for data with autocorrelation
    Estimator for the standard deviation of the mean'''
    x = np.copy(data)
    b_max = np.log2(data.size)
    b = np.empty(b_max)
    for i in range(b_max):
        b[i] = std(x) / np.sqrt(x.size)
        x = halve(x)
    return b 


def plotblocking(data):
    '''plots blocking(data) with matplotlib'''
    b = blocking(data)
    plt.plot(b)
    plt.xlabel("Block size (log2)")
    plt.ylabel("Standard deviation (naïve estimate)")
    plt.show()
    return


def binder(data):
    '''Calculates Binder Cumulant B(x) = <x^4>/(3<x^2>^2)'''
    x4 = data**4
    x2 = data**2
    d = x4.mean()
    d /= (3 * (x2.mean()**2))
    return d


def bootstrap_var(data, blocksize, n_samples):
    '''Calculates the statistical error (std) on the variance of
    a 1D array with the Bootstrap algorithm.'''
    x = np.copy(data)
    ndata = x.size
    ndiscard = ndata % blocksize
    if ndiscard > 0:
        x = x[ndiscard:]
    # In how many blocks of size blocksize can we divide our array?
    nblocks = int(ndata // blocksize)
    # Let us create an array of blocks we can choose from
    blocks = x.reshape(nblocks, blocksize)
    # This is a matrix. Each row is a sample, each element is a block (index)
    # !!! THIS OPERATION IS EXTREMELY RAM INTENSIVE !!!
    choices = np.random.randint(0, nblocks, size=(n_samples, nblocks))
    # Now let us replace each index with the block it refers to
    # And let us concatenate all blocks of each sample
    samples = blocks[choices].reshape(n_samples, nblocks*blocksize)
    # Now, axis 0 picks the sample while axis 1 goes along each sample
    # We can calculate the variance
    f_samples = samples.var(axis=1, ddof=1)
    f_std = f_samples.std(ddof=1)
    return f_std

@njit
def bootstrap_var_lessram(data, blocksize, n_samples):
    '''Calculates the statistical error (std) on the variance of
    a 1D array with the Bootstrap algorithm.'''
    x = np.copy(data)
    ndata = x.size
    ndiscard = ndata % blocksize
    if ndiscard > 0:
        x = x[ndiscard:]
    # In how many blocks of size blocksize can we divide our array?
    nblocks = int(ndata // blocksize)
    # Let us create an array of blocks we can choose from
    blocks = x.reshape(nblocks, blocksize)
    f_samples = np.empty(n_samples)
    for i in range(n_samples):
        choices = np.random.randint(0, nblocks, size=nblocks)
        sample = blocks[choices].flatten()
        f_samples[i] = var(sample)
    # Now, axis 0 picks the sample while axis 1 goes along each sample
    # We can calculate the variance
    f_std = std(f_samples)
    return f_std


def bootstrap_var_blocking(data, n_samples):
    # reblocking analysis for data with autocorrelation
    x = np.copy(data)
    b_max = int(np.log2(x.size))
    b = [bootstrap_var_lessram(x, 2**i, n_samples) for i in tqdm(range(1, b_max))]
    return np.asarray(b)

def autocorr(data):
    # returns statistical autocorrelation function of a 1D array
    x = data.copy()
    x = x - x.mean()
    aut = np.correlate(x, x, mode='full')
    aut = aut[aut.size//2:] / aut.max()
    return aut


def tint(data):
    # integrated autocorrelation time
    # returns an array, y[i] = tint integrated from 0 to i
    return np.cumsum(autocorr(data))


def thermal(f, data, xmax, step=10, **kwargs):
    '''Returns thermalization analysis of f(x)'''
    params = kwargs.get('params', None)
    measure = []
    if params is None:
        for i in range(0, xmax // step):
            measure.append(f(data[step*i:]))
    else:
        for i in range(0, xmax // step):
            measure.append(f(data[step*i:], *params))
    return np.asarray(measure)


def rt1(x):
    # round to 1st significant digit
    return round(x, -int(np.floor(np.log10(abs(x)))))


def rte(x, dx):
    # round error to measurement's 1st significant digit
    return round(x, -int(np.floor(np.log10(abs(dx)))))


def r_measurement(x, dx):
    # rounds measurement to error, assuming that error should be rounded to 1st significant digit
    # do not use indiscriminately
    return (rte(x,dx), rt1(dx))
