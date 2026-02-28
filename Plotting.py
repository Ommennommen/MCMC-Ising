import numpy as np
import matplotlib.pyplot as plt

temp = np.loadtxt('ising_temps.csv', delimiter=',')
E = np.loadtxt('ising_energy.csv', delimiter=',')
M = np.loadtxt('ising_mag.csv', delimiter=',')
C = np.loadtxt('ising_specHeat.csv', delimiter=',')
X = np.loadtxt('ising_susc.csv', delimiter=',')

fig, axs = plt.subplots(2, 2, figsize=(12, 8))

axs[0,0].scatter(temp, E, color='IndianRed')
axs[0,0].set_xlabel('Temperature'); axs[0,0].set_ylabel('Energy')

axs[0,1].scatter(temp, abs(M), color='RoyalBlue')
axs[0,1].set_xlabel('Temperature'); axs[0,1].set_ylabel('Magnetisation')

axs[1,0].scatter(temp, C, color='IndianRed')
axs[1,0].set_xlabel('Temperature'); axs[1,0].set_ylabel('Specific Heat')

axs[1,1].scatter(temp, X, color='RoyalBlue')
axs[1,1].set_xlabel('Temperature'); axs[1,1].set_ylabel('Susceptibility')

plt.tight_layout()
plt.savefig('observables.png')
plt.show()