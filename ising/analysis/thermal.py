#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 24 20:05:52 2020

Thermalization analysis of the four main observables of the Ising 2D model:
Array:for each observable O, therm_O(i) = <O> obtained discarding the first i*step measures.  
"""
import matplotlib.pyplot as plt
import numpy as np
from analysis import thermal
from tqdm import tqdm
import os
import sys

data_path = '../out/'

def chiC(data, N):
    return (N**2) * data.var(ddof=1)


# Key assumption: betas are the same for every N
# If this is not the case, the code should be adapted
Ns = np.array([20,30,40,50,60])
betas = np.array([float(strbeta[6:]) for strbeta in os.listdir(f'{data_path}20')])
betas.sort()

N = int(sys.argv[1])
beta = sys.argv[2]

e, m = np.loadtxt(f'{data_path}{N}/ising_{beta}', unpack=True)
m = np.abs(m)

i_max = 10000
step = 1000
 
therm_e = thermal(np.mean, e, i_max, step)
therm_m = thermal(np.mean, m, i_max, step)
#therm_C = thermal(chiC, e, i_max, step, params=[N])
#therm_chi = thermal(chiC, m, i_max, step, params=[N])


print(f'E: {therm_e}')
print(f'M: {therm_m}')

