import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import sys

'''
Neta = 10

N = [200, 300, 400, 500]
e1 = [18.68, 19.73, 20.09, 20.13]
de1 = [0.22, 0.17, 0.14, 0.10]

Neta = 25

N = [600, 700, 850, 1000, 1200, 1400, 1600]
e1 = [19.04, 19.31, 19.73, 19.78, 19.89, 19.88, 19.89]
de1 = [0.18, 0.08, 0.14, 0.09, 0.07, 0.06, 0.04]
'''
Neta = 5
file = sys.argv[1]

N, e, de = np.loadtxt(f'fit_out/{file}', unpack=True)


def f(x,a,b):
    return a + b*x**2 



eta = (Neta/N)
#mask = np.logical_and(N > 700, N < 1600) #N > 1000 
mask = N > 700


initp = [10,1]

popt, pcov = curve_fit(f, eta[mask], e[mask], sigma=de[mask], p0=initp, absolute_sigma=True)


chisq = (((e[mask] - f(eta[mask], *popt))/de[mask])**2).sum() / (e[mask].size - 2)

print(f'{popt[0]} +- {np.sqrt(pcov[0,0])}, χ2/ndof = {chisq}')

print(f'Valore atteso: E1 = {2*np.pi**2}')

print('Altri parametri:')

print(f'b = {popt[1]} +- {np.sqrt(pcov[1,1])}')
#print(f'c = {popt1[2]} +- {np.sqrt(pcov1[2,2])}')


#plt.xscale('log')
#plt.yscale('log')


xmin = eta[mask].min()*0.8
xmax = eta[mask].max()*1.2
ee = np.linspace(xmin, xmax, 1000)

plt.figure(1)
plt.rc('font', size=15)
plt.subplot2grid((4,1),(0,0), rowspan=3)
plt.tick_params(labelbottom=False)
plt.title(f'$\\Delta E_1$: Limite al continuo')
plt.ylabel("$\\Delta E_1$")
plt.xlim(xmin, xmax)
plt.errorbar(eta[mask], e[mask],de[mask], linestyle='', marker='o')
plt.plot(ee, f(ee, *popt))

# Residui
plt.subplot2grid((4,1),(3,0))
plt.xlim(xmin,xmax)
plt.ylim(-2.5, 2.5)
plt.ylabel('Residui')
plt.xlabel('$\\eta$')
plt.plot(eta[mask], (e[mask] - f(eta[mask], *popt))/de[mask], marker='o', linestyle='')
plt.subplots_adjust(hspace=0.14) # do not touch
plt.tight_layout()
plt.savefig(f'energy_continuum_{Neta}.pdf')
plt.show()

