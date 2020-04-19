import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
import os
import linefit_library as ls

Ls = [20, 60]


betas = []
chis = []
dchis = []
Cs = []
dCs = []



# valori teorici
n = 1
b = 1/8
g = 7/4
a = 0
bc = 0.4406868




for L in Ls:
    bb, C, dC, chi, dchi = np.loadtxt('bootstrap_out/%s' % L, unpack=True)
    vol = L**2
    betas.append(bb)
    chis.append(chi*vol)
    dchis.append(dchi*vol)
    Cs.append(C*vol)
    dCs.append(dC*vol)

# torniamo a numpy
betas = np.asarray(betas)
chis = np.asarray(chis)
dchis = np.asarray(dchis)
Cs = np.asarray(Cs)
dCs = np.asarray(dCs)


# bellurie
plt.rc('font',size=13)
plt.tight_layout()
plt.minorticks_on()



popts = []
pcovs = []

#bpc = bc - betas[chis.argmax(axis=1)]
#dbpc = 0.0025



### FIT
chi = chis[-1]
dchi = dchis[-1]
beta = betas[-1]
mask = np.logical_and(((bc - beta) > 0.03), ((bc - beta) < 0.08))
logchi = (np.log(chi))[mask]
dlogchi = (dchi/chi)[mask]
logbbc = (np.log(bc - beta))[mask]

pars, dpars, covm, chisquare = ls.linear(logbbc, logchi, dlogchi)

def chif(x, a, b):
    return x**a * np.exp(b)

xmin = 0.003
xmax = 0.1
plt.figure(1)
plt.subplot2grid((4, 1), (0, 0), rowspan=3)
plt.tick_params(labelbottom=False)
plt.title('Suscettività magnetica ($\\beta < \\beta_c$)')
plt.ylabel('$\\chi$')
plt.xscale('log')
plt.yscale('log')
plt.xlim(xmin, xmax)
plt.ylim(chis[betas < bc].min()*0.8, chis[betas < bc].max()*1.2)

for beta, chi, dchi, L in zip(betas, chis, dchis, Ls):
    plt.errorbar(bc - beta, chi, dchi, label='L = %s' % L, linestyle=' ', capsize=3)
    
bb = np.linspace(0.003,bc - betas[-1].min() , 1000)
plt.plot(bb, chif(bb, *pars), label='fit L = 60')
plt.legend()


plt.subplot2grid((4,1),(3,0))
plt.ylabel('Residui')
plt.xlabel('$\\beta_c - \\beta$')
plt.xscale('log')
plt.yscale('linear')
plt.xlim(xmin, xmax)
plt.ylim(-1,1)

chi = chis[-1]
dchi = dchis[-1]
beta = betas[-1]
plt.plot(bc - beta[mask], (np.log(chi[mask]) - np.log(chif(bc - beta[mask], *pars)))*chi[mask]/dchi[mask], linestyle=' ', marker='o')

plt.subplots_adjust(hspace=0.14) # do not touch

plt.savefig('chi_fit_rev.pdf')


########## CALORE SPECIFICO




    
### FIT
C = Cs[-1]
dC = dCs[-1]
beta = betas[-1]
mask = np.logical_and(((bc - beta) > 0.008), ((bc - beta) < 0.6))

C = C[mask]
dC = dC[mask]
logbbc = (np.log(bc - beta))[mask]

pars, dpars, covm, chisquare = ls.linear(logbbc, C, dC)


def Cf(x, a, b):
    return a*np.log(x) * b

xmin = 0.003
xmax = 0.1
plt.figure(2)
plt.subplot2grid((4, 1), (0, 0), rowspan=3)
plt.tick_params(labelbottom=False)
plt.title('Capacità termica ($\\beta < \\beta_c$)')
plt.ylabel('C')
plt.xscale('log')
plt.yscale('linear')
plt.xlim(xmin, xmax)
plt.ylim(Cs[betas < bc].min()*0.8, Cs[betas < bc].max()*1.2)

for beta, C, dC, L in zip(betas, Cs, dCs, Ls):
    plt.errorbar(bc - beta, C, dC, label='L = %s' % L, linestyle=' ', capsize=3)

bb = np.linspace(0.003,bc - betas[-1].min() , 1000)
plt.plot(bb, pars[0]*np.log(bb) + pars[1], label='fit L = 60')

plt.legend()

plt.subplot2grid((4,1),(3,0))
plt.ylabel('Residui')
plt.xlabel('$\\beta_c - \\beta$')
plt.xscale('log')
plt.yscale('linear')
plt.xlim(xmin, xmax)
plt.ylim(-1,1)

plt.plot(bc - beta[mask], (C[mask] - pars[0]*np.log(bc - beta[mask]) - pars[1])/dC[mask], linestyle=' ', marker='o')

plt.subplots_adjust(hspace=0.14) # do not touch

plt.savefig('C_fit_rev.pdf')

plt.show()
