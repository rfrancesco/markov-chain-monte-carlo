#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 24 20:05:52 2020

Correlation function analysis for estimating energy gaps
Assumption: C ~ e^{-tE} for small t -> only if <0|O|n> ~ δ_{01}
If this is not verified, an effective mass plot should be generated and studied with other tools

- Bootstrap algorithm for estimating errors

"""
import matplotlib.pyplot as plt
import numpy as np
from analysis import blocking, varblocking, std
from tqdm import tqdm
import os
import sys
from scipy.optimize import curve_fit

def eprint(x):
    print(x, file=sys.stderr)

def f(x,a,b):
    return a * np.exp(-b*x)

if len(sys.argv) == 1:
    print(f'Usage: python3 {sys.argv[0]} Neta N')
    exit(1)

Neta = int(sys.argv[1])
N = int(sys.argv[2])
folder = "tailor"

eprint(f'Reading out_{folder}/{Neta}/{N}')
data = np.loadtxt(f'out_{folder}/{Neta}/{N}.gz')

n_therm = 2*data.shape[0]//10

eprint(f'Read out_{folder}/{Neta}/{N}')
eprint(f'Discarding {n_therm} measures...')
correlators = data[n_therm:,:]

def bootstrap2d(data, blocksize, n_samples):
    x = np.copy(data)
    nmeasures = x.shape[0]
    nobservables = x.shape[1]
    eta = Neta/N
    tau = np.arange(nobservables)*eta
    mask = tau < 0.1
    ndiscard = nmeasures % blocksize
    if ndiscard > 0:
        x = x[ndiscard:]
    nblocks = int(nmeasures // blocksize)
    blocks = x.reshape(nblocks, blocksize, nobservables) 
    f_samples = np.empty(n_samples)
    for i in tqdm(range(n_samples)):
        choices = np.random.randint(0, nblocks, size=nblocks)
        sample = blocks[choices].reshape(nblocks*blocksize, nobservables)
        tempcorrelator = sample.mean(axis=0)
        popt, pcov = curve_fit(f, tau[mask], tempcorrelator[mask], p0=[1,1])
        f_samples[i] = popt[1]
    f_avg = f_samples.mean()
    f_std = std(f_samples)
    return (f_avg, f_std)

e, de = bootstrap2d(correlators, 1000, 100)

print(f'{N}\t{e}\t{de}')

