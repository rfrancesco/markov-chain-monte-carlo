import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm

Ns = [20, 60]
betas = []
chis = []
dchis = []
Cs = []
dCs = []




for N in tqdm(Ns):
    bb, C, dC, chi, dchi = np.loadtxt(f'bootstrap_out/{N}', unpack=True)
    betas.append(bb)
    vol = N**2
    chis.append(chi*vol)
    dchis.append(dchi*vol)
    Cs.append(C*vol)
    dCs.append(dC*vol)

# bellurie 
plt.rc('font',size=16)
plt.minorticks_on()

plt.figure(1)
plt.title('Suscettività magnetica')
plt.xlabel('$\\beta$')
plt.ylabel('$\\chi$')
for beta, chi, dchi, N in zip(betas, chis, dchis, Ns):
    plt.errorbar(beta, chi, dchi, label='N = %s' % N, linestyle='')
plt.legend()
plt.savefig('chi_plot.pdf')

plt.figure(2)
plt.title('Capacità termica')
plt.xlabel('$\\beta$')
plt.ylabel('C')
for beta, C, dC, N in zip(betas, Cs, dCs, Ns):
    plt.errorbar(beta, C, dC, label='N = %s' % N, linestyle='')
plt.legend()
plt.savefig('C_plot.pdf')



# valori teorici
nu = 1
b = 1/8
g = 7/4
a = 0
bc = 0.4406868



plt.figure(3)
plt.title('Magnetic susceptibility: scaling')
plt.xlabel('$(\\beta-\\beta_c)*N^{1/ \\nu}$')
plt.ylabel('$\\chi/N^{\\gamma/ \\nu }$')

for beta, chi, dchi, N in zip(betas, chis, dchis, Ns):
    dp = chi/(N**(g/nu))
    dc = dchi/(N**(g/nu))
    bb = (beta - bc)*(N**(1/nu))
    plt.errorbar(bb, dp, dc, label='N = %d' % N, marker='o', linestyle='')

plt.legend()
plt.savefig('chi_scaling.pdf')

plt.figure(4)
plt.title('Capacità termica: scaling')
plt.xlabel('$(\\beta-\\beta_c)*N^{1/\\nu}$')
plt.ylabel('$C/N^{\\alpha/\\nu}$')
for beta, C, dC, N in zip(betas, Cs, dCs, Ns):
    dp = (C - C.max())/(N**(a/nu))
    dc = (dC + dC[C.argmax()])/(N**(a/nu))
    bb = (beta - bc)*(N**(1/nu))
    plt.errorbar(bb, dp, dc, label='N = %s' % N, marker='o', linestyle='')
plt.legend()
plt.savefig('C_scaling.pdf')

plt.show()
