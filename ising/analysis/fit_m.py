import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
import os
import linefit_library as ls

Ls = [20,30,40,50,60]

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




for L in tqdm(Ls):
    bb, e, de, m, dm = np.loadtxt('blocking_out/%s' % L, unpack=True)
    betas.append(bb)
    ms.append(np.asarray(m))
    dms.append(dm)
    es.append(e)
    des.append(de)

mask = np.logical_and(((betas[-1] - bc) > 0.008), ((betas[-1] - bc) < 0.035))
# torniamo a numpy
betas = np.asarray(betas)
ms = np.asarray(ms)
dms = np.asarray(dms)
es = np.asarray(es)
des = np.asarray(des)


# bellurie 
plt.rc('font',size=13)
plt.tight_layout()
plt.minorticks_on()



popts = []
pcovs = []

# bpc = bc - betas[ms.argmax(axis=1)]
# dbpc = 0.0025



### FIT
m = ms[-1]
dm = dms[-1]
beta = betas[-1]

logm = (np.log(m))[mask]
dlogm = (dm/m)[mask]
logbbc = (np.log(beta - bc))[mask]

pars, dpars, covm, msquare = ls.linear(logbbc, logm, dlogm)

def mf(x, a, b):
    return x**a * np.exp(b)

xmin = 0.003
xmax = 0.1
plt.figure(1)
plt.subplot2grid((4, 1), (0, 0), rowspan=3)
plt.tick_params(labelbottom=False)
plt.title('Magnetizzazione per L = 60')
#plt.xlabel('$\\beta - \\beta_c$')
plt.ylabel('M')
plt.xscale('log')
plt.yscale('log')
plt.xlim(xmin, xmax)
plt.ylim(0.8*ms[betas > bc].min(), 1.2*ms[betas >bc].max())

for beta, m, dm, L in zip(betas, ms, dms, Ls):
    plt.errorbar(beta - bc, m, dm, label='L = %s' % L, linestyle=' ', capsize=3)
    
bb = np.linspace(0.003,betas[-1].max() - bc, 1000)
plt.plot(bb, mf(bb, *pars), label='fit L = 60')
plt.legend()


plt.subplot2grid((4,1),(3,0))
plt.ylabel('Residui')
plt.xlabel('$\\beta - \\beta_c$')
plt.xscale('log')
plt.yscale('linear')
plt.xlim(xmin, xmax)
#plt.ylim(-1,1)

plt.plot(betas[-1][mask] - bc, (np.log(ms[-1][mask]) - np.log(mf(betas[-1][mask] - bc, *pars)))*ms[-1][mask]/dms[-1][mask], linestyle=' ', marker='o')

plt.subplots_adjust(hspace=0.14) # do not touch

plt.savefig('m_fit.pdf')

'''
########## CALORE SPECIFICO



mask = ((betas - bc) > 0.008) & ((betas - bc) < 0.06)

    
### FIT
C = Cs[-1]
dC = dCs[-1]

C = C[mask]
dC = dC[mask]
logbbc = (np.log(betas - bc))[mask]

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

for C, dC, L in zip(Cs, dCs, Ls):
    plt.errorbar(betas - bc, C, dC, label='L = %s' % L, linestyle=' ', capsize=3)

bb = np.linspace(0.003,betas.max() - bc, 1000)
plt.plot(bb, pars[0]*np.log(bb) + pars[1], label='fit L = 60')

plt.legend()

plt.subplot2grid((4,1),(3,0))
plt.ylabel('Residui')
plt.xlabel('$\\beta - \\beta_c$')
plt.xscale('log')
plt.yscale('linear')
plt.xlim(xmin, xmax)
plt.ylim(-1,1)

plt.plot(betas[mask] - bc, (Cs[-1][mask] - pars[0]*np.log(betas[mask]-bc) - pars[1])/dCs[-1][mask], linestyle=' ', marker='o')

plt.subplots_adjust(hspace=0.14) # do not touch

plt.savefig('C_fit.pdf')
'''
plt.show()
