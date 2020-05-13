import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import os


def slab_fit(Neta, N, i):
    X, chi, dchi = np.loadtxt(f"blocking_out_slab/{Neta}/{N}", unpack=True)

    def f(x, a):
        return a * (1 - x) * x

    initpars = [1]

    popt, pcov = curve_fit(f, X, chi, sigma=dchi, p0=initpars, absolute_sigma=True)
    print(f"N = {N}")
    print(f"chi = {popt[0]} +- {np.sqrt(pcov[0,0])}")

    chisq = (((chi - f(X, *popt))/dchi)**2).sum()
    print(f'χ2 / ndof = {chisq} / 4 = {chisq/4}')

    xmin = 0
    xmax = 0.6
    ymin = 0
    ymax = 0.3
    ee = np.linspace(xmin, xmax, 1000)

    plt.rc('font',size=15)

    plt.figure(i)
    plt.subplot2grid((4,1),(0,0), rowspan=3)
    plt.xlim(xmin,xmax)
    plt.ylim(ymin,ymax)
    plt.tick_params(labelbottom=False)
    plt.title(f'Slab method: $N\\eta = {Neta}, N = {N}$')
    plt.ylabel("$\\chi_{x,\\eta}$")
    plt.errorbar(X, chi, dchi, linestyle='', marker='o')
    plt.plot(ee, f(ee, *popt), label='fit')

    # Residui
    plt.subplot2grid((4,1),(3,0))
    plt.xlim(xmin,xmax)
    plt.ylim(-2,2)
    plt.ylabel('Residui')
    plt.xlabel('$x$')
    plt.plot(X, (chi - f(X, *popt))/dchi, marker='o', linestyle='')
    plt.subplots_adjust(hspace=0.14) # do not touch
    plt.tight_layout()
    plt.savefig(f'slab_fit_{Neta}_{N}.pdf')
    #plt.show()
    return (popt[0], np.sqrt(pcov[0,0]))


if __name__ == "__main__":
    Neta = 5
    Ns = np.array([1050, 1150, 1250, 1350, 1500])
    eta = Neta / Ns
    chi = []
    dchi = []
    for i, N in enumerate(Ns):
        c, dc = slab_fit(Neta, N, i)
        chi.append(c)
        dchi.append(dc)


    def g(x, a, b):
        return a + b*x**2

    initpars = [1, 1]
    popt, pcov = curve_fit(g, eta, chi, sigma=dchi, p0=initpars, absolute_sigma=True)
    print(f'{popt[0]} +- {np.sqrt(pcov[0,0])}')

    chisq = (((chi - g(eta, *popt))/dchi)**2).sum()
    print(f'χ2 / ndof = {chisq} / 4 = {chisq/4}')

    xmin = 0.002
    xmax = eta.max()*1.2
    ee = np.linspace(xmin, xmax, 1000)

    plt.figure(Ns.size + 1)
    plt.rc('font', size=15)
    plt.subplot2grid((4, 1), (0, 0), rowspan=3)
    plt.title(f"Slab method: $N\\eta = {Neta}$, Limite $\\eta \\to 0$")
    plt.xlim(xmin,xmax)
    plt.tick_params(labelbottom=False)
    plt.ylabel("$\\chi_\\eta$")
    plt.xlabel("")
    plt.errorbar(eta, chi, dchi, linestyle='', marker='o')
    plt.plot(ee, g(ee, *popt), label='fit')

    # Residui
    plt.subplot2grid((4, 1), (3, 0))
    plt.xlim(xmin, xmax)
    plt.ylim(-2, 2)
    plt.ylabel('Residui')
    plt.xlabel('$\\eta$')
    plt.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
    plt.plot(eta, (chi - g(eta, *popt))/dchi, marker='o', linestyle='')
    plt.subplots_adjust(hspace=0.14) # do not touch
    plt.tight_layout()
    plt.savefig(f'slab_fit_{Neta}_continuum.pdf')
    plt.show()
