import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
import os
import linefit_library as ls

Ls = [20,30, 40, 50, 60]
Ls = np.array(Ls)

betas = []
ms = []
dms = []
es = []
des = []



# valori teorici
n = 1
b = 1/8
g = 7/4
a = 0
bc = 0.4406868
beta = bc

N, C, dC, chi, dchi = np.loadtxt('bootstrap_out', unpack=True)
C = C * N**2
dC = dC * N**2
chi = chi * N**2
dchi = dchi * N**2


# bellurie
plt.rc('font',size=13)
plt.tight_layout()
plt.minorticks_on()



popts = []
pcovs = []

#bpc = bc - betas[chis.argmax(axis=1)]
#dbpc = 0.0025



### FIT
logchi = (np.log(chi))
dlogchi = (dchi/chi)
logN = (np.log(N))

pars, dpars, covm, chisquare = ls.linear(logN, logchi, dlogchi)

def chif(x, a, b):
    return x**a * np.exp(b)

#xmin = 0.003
#xmax = 0.1
plt.figure(1)
plt.subplot2grid((4, 1), (0, 0), rowspan=3)
plt.tick_params(labelbottom=False)
plt.title('Suscettività magnetica a $\\beta_c$')
plt.xlabel('$N')
plt.ylabel('$\\chi$')
plt.xscale('log')
plt.yscale('log')
#plt.xlim(xmin, xmax)

plt.errorbar(N, chi, dchi, label='', linestyle=' ', capsize=3)
    
bb = np.linspace(0.8*N.min(),1.2*N.max(), 1000)
plt.plot(N, chif(N, *pars), label='fit')
plt.legend()


plt.subplot2grid((4,1),(3,0))
plt.ylabel('Residui')
plt.xlabel('$N$')
plt.xscale('log')
plt.yscale('linear')
#plt.xlim(xmin, xmax)
plt.ylim(-1,1)
plt.plot(N, (np.log(chi) - np.log(chif(N, *pars)))*chi/dchi, linestyle=' ', marker='o')
plt.subplots_adjust(hspace=0.14) # do not touch

plt.savefig('chi_fs_fit.pdf')
plt.show()
'''
########## CALORE SPECIFICO




    
### FIT
C = Cs[-1]
dC = dCs[-1]
beta = betas[-1]
mask = np.logical_and(((beta - bc) > 0.008), ((beta - bc) < 0.6))

C = C[mask]
dC = dC[mask]
logbbc = (np.log(beta - bc))[mask]

pars, dpars, covm, chisquare = ls.linear(logbbc, C, dC)


def Cf(x, a, b):
    return a*np.log(x) * b

xmin = 0.003
xmax = 0.1
plt.figure(2)
plt.subplot2grid((4, 1), (0, 0), rowspan=3)
plt.tick_params(labelbottom=False)
plt.title('Capacità termica per L = 60')
plt.ylabel('C')
plt.xscale('log')
plt.yscale('linear')
plt.xlim(xmin, xmax)

for beta, C, dC, L in zip(betas, Cs, dCs, Ls):
    plt.errorbar(beta - bc, C, dC, label='L = %s' % L, linestyle=' ', capsize=3)

bb = np.linspace(0.003,betas[-1].max() - bc, 1000)
plt.plot(bb, pars[0]*np.log(bb) + pars[1], label='fit L = 60')

plt.legend()

plt.subplot2grid((4,1),(3,0))
plt.ylabel('Residui')
plt.xlabel('$\\beta - \\beta_c$')
plt.xscale('log')
plt.yscale('linear')
plt.xlim(xmin, xmax)
plt.ylim(-1,1)

plt.plot(beta[mask] - bc, (C[mask] - pars[0]*np.log(beta[mask]-bc) - pars[1])/dC[mask], linestyle=' ', marker='o')

plt.subplots_adjust(hspace=0.14) # do not touch

plt.savefig('C_fit.pdf')

plt.show()
'''
