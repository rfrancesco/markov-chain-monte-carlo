import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
import os

data_path = '../out/'

Ns = np.array([20,30, 40, 50, 60])

betas = []
E = []
dE = []
M = []
dM = []

# valori teorici
nu = 1
b = 1/8
g = 7/4
a = 0
bc = 0.4406868


for N in tqdm(Ns):
    bb, ee, dee, mm, dmm = np.loadtxt('blocking_out/%s' % N, unpack=True)
    betas.append(bb)
    E.append(ee)
    dE.append(dee)
    M.append(mm)
    dM.append(dmm)

plt.figure(1)
plt.title('Densità di energia')
plt.xlabel('$\\beta$')
plt.ylabel('$E/N^2$')
for i, N in enumerate(Ns):
    plt.errorbar(betas[i], E[i], dE[i], label='N = %s' % N, linestyle=' ', marker='+')
plt.legend()
plt.savefig('e_plot.pdf')

plt.figure(2)
plt.title('Magnetizzazione')
plt.xlabel('$\\beta$')
plt.ylabel('$M/N^2$')
for i, N in enumerate(Ns):
    plt.errorbar(betas[i], M[i], dM[i], label='N = %s' % N, linestyle=' ', marker='+')
plt.legend()
plt.savefig('m_plot.pdf')

plt.figure(3)
plt.title('Magnetizzazione: scaling')
plt.xlabel('$\\beta$')
plt.ylabel('$M/N^2$')
for i, N in enumerate(Ns):
    bb = (betas[i] - bc)*N**(1/nu)
    mm = M[i]/N**(-b/nu)
    dmm = dM[i]/N**(-b/nu)
    plt.errorbar(bb, mm, dmm, label='N = %s' % N, linestyle=' ', marker='+')
plt.legend()
plt.savefig('m_scaling.pdf')
plt.show()
