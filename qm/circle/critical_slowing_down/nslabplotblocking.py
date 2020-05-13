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
    print(f'Usage: python3 {sys.argv[0]} Neta N')
    exit(1)

folder = "slab"
    
Neta = int(sys.argv[1])
N = int(sys.argv[2])

data = np.loadtxt(f'blocking_raw_{folder}/{Neta}/{N}')

Ns = data[:,0]
Xs = data[:,1]
chis = data[:,2]
blocks = data[:,3:]


for i, (N, X, chi, block) in enumerate(zip(Ns, Xs, chis, blocks)):
    chiblock = block[:block.size//2]
    taublock = block[block.size//2:]
    plt.subplot(2,1,1)
    plt.title(f'Neta = {Neta}, N = {N}, x = {X} => $\\eta$ = {Neta/N:.3f},  $\\chi$ = {chi}')
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


