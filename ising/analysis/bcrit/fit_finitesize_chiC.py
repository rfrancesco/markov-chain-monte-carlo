import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
import os
import linefit_library as ls

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
plt.rc('font',size=15)
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

xmin = N.min() * 0.9
xmax = N.max() * 1.1
plt.figure(1)
plt.subplot2grid((4, 1), (0, 0), rowspan=3)
plt.tick_params(labelbottom=False)
plt.title('Suscettività magnetica a $\\beta_c$')
plt.xlabel('$N')
plt.ylabel('$\\chi$')
plt.xscale('log')
plt.yscale('log')
plt.xlim(xmin, xmax)

plt.errorbar(N, chi, dchi, label='', linestyle=' ', capsize=3)
    
bb = np.linspace(0.8*N.min(),1.2*N.max(), 1000)
plt.plot(N, chif(N, *pars), label='fit')
plt.legend()


plt.subplot2grid((4,1),(3,0))
plt.ylabel('Residui')
plt.xlabel('$N$')
plt.xscale('log')
plt.yscale('linear')
plt.xlim(xmin, xmax)
plt.ylim(-2,2)
plt.plot(N, (np.log(chi) - np.log(chif(N, *pars)))*chi/dchi, linestyle=' ', marker='o')
plt.tight_layout()
plt.subplots_adjust(hspace=0.1) # do not touch

plt.savefig('../../relazione/figure/chi_fs_fit.pdf')

########## CALORE SPECIFICO




    
### FIT
C = C
dC = dC
logN = (np.log(N))

pars, dpars, covm, chisquare = ls.linear(logN, C, dC)


def Cf(x, a, b):
    return a*np.log(x) * b

plt.figure(2)
plt.subplot2grid((4, 1), (0, 0), rowspan=3)
plt.tick_params(labelbottom=False)
plt.title('Calore specifico a $\\beta_c$')
plt.ylabel('C')
plt.xscale('log')
plt.yscale('linear')
plt.xlim(xmin, xmax)

plt.errorbar(N, C, dC, linestyle=' ', capsize=3)

bb = np.linspace(0.8*N.min(),1.2*N.max(), 1000)
plt.plot(bb, pars[0]*np.log(bb) + pars[1], label='fit')

plt.legend()

plt.subplot2grid((4,1),(3,0))
plt.ylabel('Residui')
plt.xlabel('$N$')
plt.xscale('log')
plt.yscale('linear')
plt.xlim(xmin, xmax)
plt.ylim(-2,2)

plt.plot(N, (C - pars[0]*np.log(N) - pars[1])/dC, linestyle=' ', marker='o')

plt.tight_layout()
plt.subplots_adjust(hspace=0.1) # do not touch

plt.savefig('../../relazione/figure/C_fs_fit.pdf')


plt.show()

