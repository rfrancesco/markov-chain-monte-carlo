import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import os


def continuum_fit(Neta, alg):
    N, chi, dchi, tau = np.loadtxt(f"blocking_out_{alg}/{Neta}", unpack=True)
    eta = Neta/N
    eta2 = eta**2

    mask = eta < 0.03

    def f(x, a, b):
        return a + b*x**2

    initpars = [1,1]

    popt, pcov = curve_fit(f, eta[mask], chi[mask], sigma=dchi[mask], p0=initpars, absolute_sigma=True)
    print(f"chi = {popt[0]} +- {np.sqrt(pcov[0,0])}")

    chisq = ((chi[mask] - f(eta[mask], *popt))/dchi[mask])**2 
    chisq = chisq.sum() / (chi[mask].size - 2)
    print(f"χ2 / ndof = {chisq}")

    xmin = 0
    xmax = eta[mask].max() * 1.2
    ymin = chi[mask].min() * 0.99
    ymax = chi[mask].max() * 1.01
    ee = np.linspace(xmin, xmax, 1000)

    plt.figure(1)
    plt.rc('font',size=15)

    plt.subplot2grid((4,1),(0,0), rowspan=3)
    plt.xlim(xmin,xmax)
    plt.ylim(ymin,ymax)
    plt.tick_params(labelbottom=False)
    plt.title(f'Misura di $\\chi$ per $N\\eta = {Neta}$')
    plt.ylabel("Misura")
    plt.errorbar(eta, chi, dchi, linestyle='')
    plt.plot(ee, f(ee, *popt), label='fit')

    # Residui
    plt.subplot2grid((4,1),(3,0))
    plt.xlim(xmin,xmax)
    plt.ylim(-2,2)
    plt.ylabel('Residui')
    plt.xlabel('$\\eta$')
    plt.plot(eta[mask], (chi[mask] - f(eta[mask], *popt))/dchi[mask], marker='o', linestyle='')
    plt.subplots_adjust(hspace=0.14) # do not touch
    plt.tight_layout()
    plt.savefig(f'continuum_fit_{alg}_{Neta}.pdf')
    plt.show()


if __name__ == "__main__":
    algs = ["tailor"]
    Netas = [5,10,15]
    for Neta in Netas:
        continuum_fit(Neta, "tailor")
