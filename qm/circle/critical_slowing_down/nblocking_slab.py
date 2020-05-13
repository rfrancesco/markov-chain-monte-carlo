#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 24 20:05:52 2020

"""
import matplotlib.pyplot as plt
import numpy as np
from analysis import blocking, varblocking
from tqdm import tqdm
import os
import sys
#from p_tqdm import p_map

if len(sys.argv) == 1:
    print(f'Usage: python3 {sys.argv[0]} Neta N')
    exit(1)

Neta = int(sys.argv[1])
N = int(sys.argv[2])

#Ns = [int(strN) for strN in os.listdir(f'out_{folder}/{Neta}/')]
Ns = [1350, 1450]
Xs = [0.1, 0.2, 0.3, 0.4, 0.5]


folder = "slab"

def getblocks(N,x):
    q = np.loadtxt(f'out_{folder}/{Neta}/{N}/{x}', unpack=True)
    print(f'Read out_{folder}/{Neta}/{N}/{x}')
    q = q[100000:]
    q2 = q**2
    chi = q2.mean() / Neta
    chiblock = blocking(q2) / Neta
    var_naive = q2.var(ddof=1) / q2.size
    tau = varblocking(q2) / var_naive 
    blocks = np.concatenate(([N], [x], [chi], chiblock, tau))
    return blocks

#results = p_map(getblocks, Ns)
results = [getblocks(N,X) for X in tqdm(Xs)]

np.savetxt(f'blocking_raw_slab/{Neta}/{N}', np.vstack(results))
