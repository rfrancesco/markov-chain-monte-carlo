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

data = np.loadtxt(f'bootstrap_raw/{N}')

betas = data[:,0]
evars = data[:,1]
mvars = data[:,2]
des = []
dms = []
blocks = data[:,3:]

for beta, evar, mvar, block in zip(betas, evars, mvars, blocks):
    eblock = block[:block.size//2]
    mblock = block[block.size//2:]
    #print(f'{beta} -> e = {evar}, m = {mvar}; \neblock = {eblock}; \nmblock = {mblock}')
    plt.subplot(2,1,1)
    plt.tight_layout()
    plt.title(f'b = {beta}, Var[E] = {evar:.7f}, Var[M] = {mvar:.7f}')
    plt.ylabel('dVar[E]')
    plt.grid(axis='y')
    plt.yticks(np.linspace(eblock.min(), eblock.max(), 10))
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.plot(eblock)
    plt.subplot(2,1,2)
    plt.tight_layout()
    plt.ylabel('dVar[M]')
    plt.grid(axis='y')
    plt.yticks(np.linspace(mblock.min(), mblock.max(), 10))
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    #plt.tight_layout()
    plt.plot(mblock)
    print(f"(N, beta) = ({N}, {beta})")
    print(f"eblock = {eblock}")
    print(f"mblock = {mblock}")
    plt.show(block=False)
    confirmed = False
    while not confirmed:
        de = inputfloat("dVarE? ")
        dm = inputfloat("dVarM? ")
        print(f"{beta}: Var[E] = {evar} +- {de}, M = {mvar} +- {dm}")
        answer = ' '
        while (answer != 'y') and (answer != 'n'):
            answer = input('Ok? (y/n)')
        if answer == 'y':
            confirmed = True
    des.append(de)
    dms.append(dm)
    print(f"Saved {beta}: Var[E] = {evar} +- {de}, Var[M] = {mvar} +- {dm}")
    plt.close()

np.savetxt(f"bootstrap_out/{N}", np.column_stack([betas, evars, des, mvars, dms]), header="beta\tVarE\tdVarE\tVarM\tdVarM\t")


    


