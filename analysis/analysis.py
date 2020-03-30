import numpy as np
import random as ran
import matplotlib.pyplot as plt
from tqdm import tqdm

def sample_var(data):
    d = data.var(ddof=1)
    return d

def sample_std(data):
    return data.std(ddof=1)

def mean_var(data):
    return (sample_var(data) / np.size(data))

def mean_std(data):
    return np.sqrt(mean_var(data))

def halve(data):
    # from 1D array x, returns array y of size x.size/2
    # where y[i] = (x[2i] + x[2i + 1]) /2
    # (a_1, ..., a_n) -> (a_1 + a_2, a_3 + a_4, ...)/2
    s = data.size
    w = np.copy(data)
    if ((s % 2) != 0):
        w = w[1:]
        s -= 1
    s /= 2
    s = int(s)
    return (w.reshape(s, 2).sum(axis=1) / 2)

        
def blocking(data):
    # reblocking analysis for data with autocorrelation
    # returns b[i] = std of data.mean() after averaging on blocks of size 2**i
    # obviously, b[0] = data.std(ddof=1)
    b = []
    x = np.copy(data)
    b.append(mean_std(x))
    for i in range(0, int(np.log2(np.size(data))) - 1):
        x = halve(x)
        b.append(mean_std(x))
    return np.asarray(b)


def plotblocking(data):
    # plots blocking(data) with matplotlib
    b = blocking(data)
    plt.plot(b)
    plt.xlabel("Block size (log2)")
    plt.ylabel("Standard deviation (naïve estimate)")
    plt.show()
    return


def binder(data):
    # binder cumulant
    x4 = data**4
    x2 = data**2
    d = x4.mean() 
    d /= (3 * (x2.mean()**2))
    return d


def optfbootstrap_wbins(f, data, blocksize, n_samples):
    # implementation of bootstrap algorithm
    # very slow code! beware!
    x = np.copy(data)
    samples = []
    ndata = x.size 
    ndiscard = ndata % blocksize
    if (ndiscard != 0):
        x = x[ndiscard:]
    nblocks = int(ndata // blocksize) 
    blocks = x.reshape(nblocks, blocksize)
    for i in tqdm(range(0, n_samples)):
        sample = []
        for j in range(0, nblocks):
            x = np.random.randint(0, nblocks)
            sample.append(blocks[x])
        samples.append(f(sample))
    f_std = np.std(samples, ddof=1)
    return f_std

def xfbootstrap_wbins(f, data, blocksize, n_samples):
    # extra fast bootstrap w numpy
    # going to replace optfbootstrap eventually
    x = np.copy(data)
    ndata = x.size 
    ndiscard = ndata % blocksize
    if (ndiscard != 0):
        x = x[ndiscard:]
    nblocks = int(ndata // blocksize) 
    blocks = x.reshape(nblocks, blocksize)
    idx = np.random.randint(0, nblocks, size=(n_samples, nblocks))
    samples = np.apply_along_axis(f, 1, blocks[idx])
    f_std = np.std(samples, ddof=1)
    return f_std


def autocorr(x):
    # returns statistical autocorrelation function of a 1D array
    tmp = x.copy()
    tmp = tmp - tmp.mean()
    aut = np.correlate(tmp, tmp, mode='full')
    aut = aut[aut.size//2:] / aut.max()
    return aut

def tint(x):
    # integrated autocorrelation time
    # returns an array, y[i] = tint integrated from 0 to i
    return np.cumsum(autocorr(x))

def model_y2(neta):
    y = 0.5
    y += 1 / (np.exp(neta) - 1)
    return y

def thermal(f, data, xmax, step=10):
    # Thermalization analysis of f(data) from 0 to xmax
    measure = []
    for i in range(0, xmax // step):
        measure.append(f(data[step*i:]))
    return np.asarray(measure)

def rt1(x): 
    #round to 1st significant digit
    return round(x, -int(np.floor(np.log10(abs(x)))))

def rte(x,dx):  
    # round error to measurement's 1st significant digit
    return round(x, -int(np.floor(np.log10(abs(dx)))))

def r_measurement(x,dx):
    # rounds measurement to error, assuming that error should be rounded to 1st significant digit
    # do not use indiscriminately
    return (rte(x,dx), rt1(dx))
