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
    print(f'Usage: python3 {sys.argv[0]} N')
    exit(1)
    
    
def inputfloat(prompt):
    while True:
        try:
            x = float(input(prompt))
            break
        except:
            print("Illegal value!")
    return x

init = 0
N = sys.argv[1]

data = np.loadtxt(f'blocking_raw/{N}')

betas = data[:,0]
es = data[:,1]
ms = data[:,2]
des = []
dms = []
blocks = data[:,3:]

for i, (beta, e, m, block) in enumerate(zip(betas, es, ms, blocks)):
    eblock = block[:block.size//2]
    mblock = block[block.size//2:]
    #print(f'{beta} -> e = {e}, m = {m}; \neblock = {eblock}; \nmblock = {mblock}')
    plt.subplot(2,1,1)
    plt.title(f'b = {beta}, E = {e:.5f}, M = {m:.5f}')
    plt.ylabel('dE')
    plt.grid(axis='y')
    plt.yticks(np.linspace(eblock.min(), eblock.max(), 10))
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.plot(eblock)
    plt.subplot(2,1,2)
    plt.ylabel('dM')
    plt.grid(axis='y')
    plt.yticks(np.linspace(mblock.min(), mblock.max(), 10))
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.plot(mblock)
    print(f"(N, beta) = ({N}, {beta})")
    print(f"eblock = {eblock}")
    print(f"mblock = {mblock}")
    plt.show(block=False)
    confirmed = False
    while not confirmed:
        de = inputfloat("dE?")
        dm = inputfloat("dM?")
        print(f"{beta}: E = {e} +- {de}, M = {m} +- {dm}")
        answer = ' '
        while (answer != 'y') and (answer != 'n'):
            answer = input('Ok? (y/n)')
        if answer == 'y':
            confirmed = True
    des.append(de)
    dms.append(dm)
    print(f"Saved {beta}: E = {e} +- {de}, M = {m} +- {dm}")
    plt.close()

np.savetxt(f"blocking_out/{N}", np.column_stack([betas, es, des, ms, dms]), header="beta\te\tde\tm\tdm\t")

