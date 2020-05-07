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
from cycler import cycler

if len(sys.argv) == 1:
    print(f'Usage: python3 {sys.argv[0]} (local/tailor) Neta')
    exit(1)

Neta = int(sys.argv[2])
folder = sys.argv[1]

Ns = [int(strN) for strN in os.listdir(f'out_{folder}/{Neta}/')]
Ns.sort()

algoritmo = {
        'local': "Metropolis locale",
        'tailor': "Metropolis non-locale"
}

# bellurie
plt.rc('font',size=13)
plt.tight_layout()


plt.figure(1)
plt.title(f'Storia Montecarlo: $N\\eta = {Neta}$, {algoritmo.get(folder, "errore")}')
plt.xlabel('Misure (1 ogni 10 spazzate)')
plt.ylabel('Q')

#plt.xlim(500000,510000)

#csd_cycler = cycler(color=['b','y','r','m','k'])
#plt.rc('axes', prop_cycle=csd_cycler)

for N in tqdm(Ns):
    q = np.loadtxt(f'out_{folder}/{Neta}/{N}')
    plt.plot(q, label=f'N = {N}', rasterized=True)

plt.legend()

plt.savefig(f'csd_plot_{folder}_{Neta}.pdf')
