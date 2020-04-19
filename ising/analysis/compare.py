#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 24 20:05:52 2020

Plots Monte Carlo history for the Ising model
For N, beta, compares hot start vs cold start
"""
import matplotlib.pyplot as plt
import numpy as np
plt.rcParams['figure.figsize'] = [10, 10]


def compare(n, beta, t_max=5000):
    '''Compares MC histories of two simulations with different initial conditions. 
    Warning: parameters must be strings for trailing zeroes in beta to be preserved.
    '''
    inits = [0, 1]

    cold, hot = [np.loadtxt(f'{init}/{n}/ising_{beta}', unpack=True) for init in inits]

    e_cold, m_cold = cold
    e_hot, m_hot = hot

    m_hot = np.abs(m_hot)
    m_cold = np.abs(m_cold)

    plt.subplot(2, 1, 1)
    plt.title(f'Monte Carlo history: L = {n}, beta = {beta}')
    plt.ylabel('E')
    plt.xlim(0, t_max)
    plt.plot(e_cold, color='blue')
    plt.plot(e_hot, color='red')

    plt.subplot(2, 1, 2)
    plt.ylabel('M')
    plt.xlim(0, t_max)
    plt.ylim(0, 1)
    plt.plot(m_cold, color='blue')
    plt.plot(m_hot, color='red')
    plt.show()

if __name__ == '__main__':
    import sys
    if len(sys.argv) < 3:
        print(f"Usage: python {sys.argv[0]} N beta (t_max) OR import")
        sys.exit(1)

    n = sys.argv[1]
    beta = sys.argv[2]
    if len(sys.argv) == 4:
        t_max = int(sys.argv[3])
        compare(n, beta, t_max)
    else:
        compare(n, beta)
    
