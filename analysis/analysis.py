import numpy as np
import random as ran
import matplotlib.pyplot as plt
from tqdm import tqdm

# Numba does not support arguments in np.var(), np.std()
# Therefore, var() and std() are equivalent to np.var(ddof=1) and np.std(ddof=1)
# (njit disabled in the code)

def var(data):
    '''Numba-friendly shorthand for np.var(data, ddof=1).
    Estimator for the variance of a data set'''
    n = data.size
    return data.var() * n/(n-1)


def mean_var(data):
    '''Variance of the mean of a data set'''
    return var(data) / data.size

def std(data):
    '''Numba-friendly shorthand for np.std(data, ddof=1).
    Estimator for the standard deviation of a data set'''
    return np.sqrt(var(data))

def mean_std(data):
    '''Standard deviation of the mean of a data set'''
    return np.sqrt(mean_var(data))


# Blocking error estimation utilities

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

def blocking(data):
    '''Blocking analysis for data with autocorrelation
    Estimator for the standard deviation of the mean'''
    x = np.copy(data)
    b_max = int(np.log2(data.size))
    b = np.empty(b_max)
    for i in range(b_max):
        b[i] = mean_std(x)
        x = halve(x)
    return b 

def varblocking(data):
    '''Blocking analysis for data with autocorrelation
    Estimator for the variance of the mean'''
    x = np.copy(data)
    b_max = int(np.log2(data.size))
    b = np.empty(b_max)
    for i in range(b_max):
        b[i] = mean_var(x)
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


# Bootstrap error estimation utilities

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
    f_samples = np.empty(n_samples)
    for i in range(n_samples):
        choices = np.random.randint(0, nblocks, size=nblocks)
        sample = blocks[choices].flatten()
        f_samples[i] = var(sample)
    # Now, axis 0 picks the sample while axis 1 goes along each sample
    # We can calculate the standard deviation
    # Note that it is NOT mean_std! Otherwise, we could reduce our error
    # by increasing n_samples
    # (See Newman, Barkema: Monte-Carlo methods in Statistical Physics)
    f_std = std(f_samples)
    return f_std


def bootstrap_var_blocking(data, n_samples):
    '''Calculates the statistical error (std) on the variance of
    a 1D array of autocorrelated data, with the Bootstrap algorithm
    with blocks of increasing size'''
    x = np.copy(data)
    b_max = int(np.log2(x.size))
    b = [bootstrap_var(x, 2**i, n_samples) for i in tqdm(range(1, b_max))]
    return np.asarray(b)

# Data autocorrelation utilities

def autocorr(data):
    '''Returns statistical autocorrelation of a 1D array'''
    x = data.copy()
    x = x - x.mean()
    aut = np.correlate(x, x, mode='full')
    aut = aut[aut.size//2:] / aut.max()
    return aut


def tint(data):
    '''Estimate of autocorrelation time. Returns integrated
    autocorrelation time tau_int, integrated from 0 to i.'''
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

# Rounding utilities

def order(x):
    '''Gets order of magnitude exponent of x (x ~ 10**order(x))'''
    return int(np.floor(np.log10(x)))

def rt1(x):
    '''Rounds x to the first significant digit'''
    return round(x, - order(x))

def rtn(x, n):
    '''Rounds x to the nth significant digit'''
    return round(x, - order(x) + (n - 1))

def rwe(x, dx, n):
    '''Rounds (x \pm dx) to nth significant digit of dx.'''
    return (round(x, -order(dx)+(n-1)), round(dx, -order(dx)+(n-1)))
