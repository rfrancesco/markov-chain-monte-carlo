import sys
import analysis as a
import numpy as np
import matplotlib.pyplot as plt

plot = True
data = np.loadtxt(sys.argv[1], skiprows=1)
thermal = int(sys.argv[2])
observable_names = ["P_t", "P_x", "P_tP_x", "P_tP_x^\\dag"]

print(f"Analyzing file {sys.argv[1]}, skipping {thermal}")

# Transpose the matrix (so that observable[i] = [measures of O_i]
# and discard {thermal} measures 
observables = data[thermal:, :].T

# Compute the module |O_i|
observables = observables**2
# This command sums over consecutive observables (rows): [A,B,C,D] -> [A+B, C+D]
observables = observables.reshape(observables.shape[0]//2, 2, observables.shape[1]).sum(axis=1)
observables = np.sqrt(observables)

# If observable_names is not set correctly, fallback to i
if (observables.shape[0] != observable_names):
    observable_names = range(0, observables.shape[0])

print(f"Number of (observables, measures): {observables.shape}")

print("## Averages")
observable_means = observables.mean(axis=1)
print(observable_means)

print("##Blocking")
for observable in observables:
    print(a.blocking(observable))

if plot:
    for (i, observable) in enumerate(observables):
        plt.figure(i)
        plt.title(f"Blocking for ${observable_names[i]}$ = {observable_means[i]:e}", alpha=0.7)
        plt.plot(a.blocking(observable))
    plt.show()
