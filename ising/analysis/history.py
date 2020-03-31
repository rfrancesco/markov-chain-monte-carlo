#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 24 20:05:52 2020

Plots Monte Carlo history for the Ising model
"""
import matplotlib.pyplot as plt
import numpy as np
from analysis import thermal
from tqdm import tqdm
import os

data_path = '../out/'

# Key assumption: betas are the same for every N
# If this is not the case, the code should be adapted
Ns = np.array([20,30,40,50,60])
betas = np.array([float(strbeta[6:]) for strbeta in os.listdir(f'{data_path}20')])
betas.sort()

for N in tqdm(Ns):
    for beta in tqdm(betas):
        e, m = np.loadtxt(f'{data_path}%d/ising_%s' % (N, beta), unpack=True)

        plt.subplot(2,1,1)
        plt.title(f'Monte Carlo history: N = {N}, beta = {beta}')
        plt.xlim(0,500)
        plt.ylabel('E')
        plt.plot(e)

        plt.subplot(2,1,2)
        plt.xlim(0,500)
        plt.ylabel('M')
        plt.plot(m)
        plt.show()
    

