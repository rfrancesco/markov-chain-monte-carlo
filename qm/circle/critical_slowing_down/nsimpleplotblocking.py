#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 24 20:05:52 2020

"""
import matplotlib.pyplot as plt
import numpy as np
import os
import sys

if len(sys.argv) == 1:
    print(f'Usage: python3 {sys.argv[0]} (local/tailor) Neta')
    exit(1)
    
    
Neta = int(sys.argv[2])
folder = sys.argv[1]

data = np.loadtxt(f'blocking_raw_{folder}/{Neta}_2')

Ns = data[:,0]
chis = data[:,1]
blocks = data[:,2:]


for i, (N, chi, block) in enumerate(zip(Ns, chis, blocks)):
    chiblock = block[:block.size//2]
    taublock = block[block.size//2:]
    plt.subplot(2,1,1)
    plt.title(f'Neta = {Neta}, N = {N} => $\\eta$ = {Neta/N},  $\\chi$ = {chi}')
    plt.ylabel('d$\\chi$')
    plt.grid(axis='y')
    plt.yticks(np.linspace(chiblock.min(), chiblock.max(), 10))
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.plot(chiblock)
    plt.subplot(2,1,2)
    plt.ylabel('$\\tau$')
    plt.grid(axis='y')
    plt.yticks(np.linspace(taublock.min(), taublock.max(), 10))
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.plot(taublock)
    plt.show()


