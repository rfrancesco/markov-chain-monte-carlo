#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 24 20:05:52 2020

"""
import matplotlib.pyplot as plt
import numpy as np
from analysis import bootstrap_var_blocking
from tqdm import tqdm
import os
import sys

if len(sys.argv) == 1:
    print(f'Usage: python3 {sys.argv[0]} N')
    exit(1)

init = 0
n_samples = 1000
N = sys.argv[1]

strbetas = [strbeta[6:] for strbeta in os.listdir(f'{init}/{N}/')]
strbetas.sort()
betas = np.array([float(strbeta) for strbeta in strbetas])


def getblocks(strbeta):
    e, m = np.loadtxt(f'{init}/{N}/ising_{strbeta}', unpack=True)
    e = e[-100000:]
    m = m[-100000:]
    print(f'Read {N}/ising_{strbeta}')
    m = np.abs(m)
    eblock = bootstrap_var_blocking(e, n_samples)
    mblock = bootstrap_var_blocking(m, n_samples)
    print(f'{strbeta} -> eblock = {eblock}, mblock = {mblock}')
    blocks = np.concatenate(([float(strbeta)], eblock, mblock))
    return blocks

np.savetxt(f'bootstrap/{N}', np.vstack([getblocks(strbeta) for strbeta in tqdm(strbetas)]))
