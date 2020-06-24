import numpy as np
import matplotlib.pyplot as plt
from itertools import product

'''Analysis of the critical slowing down
   Plotting (1/Œ , t_int). '''

Netas = [5, 10, 15]
algs = ["local", "tailor"]

plt.rc('font',size=13)
plt.tight_layout()
plt.figure(1)
plt.yscale('log')
plt.title('Confronto di $\\tau_\\chi$ tra algoritmo locale e non-locale')
plt.xlabel('$1/\\eta$')
plt.ylabel('$\\tau_\\chi$ [spazzate]')

markers = {
    "local" : "o",
    "tailor" : "D"
}

for (alg, Neta) in product(algs, Netas):
    print(f'{Neta}, {alg}')
    N, chi, dchi, tau = np.loadtxt(f"blocking_out_{alg}/{Neta}", unpack=True)
    tau = tau*10 / 2 # fattore /2: mi sono dimenticato di dividere, quindi in realt√† sono tutti sballati di un fattore 2
    plt.plot(N/Neta, tau, linestyle='', marker=markers[alg], label=f'$N\\eta$ = {Neta}, {alg}')

plt.legend()
plt.savefig("csd_autocorr.pdf")
plt.savefig("png_csd_autocorr.png")
plt.show()
