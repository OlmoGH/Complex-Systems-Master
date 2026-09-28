import numpy as np
import matplotlib.pyplot as plt
import os

directorio = os.path.dirname(os.path.abspath(__file__))

data_energy = np.loadtxt(directorio + "/out/Energy.txt")
energy = data_energy[:-1]
c_H = data_energy[-1]

data_magnetization = np.loadtxt(directorio + "/out/Magnetization.txt")
magnetization = data_magnetization[:-1]
chi = data_magnetization[-1]

plt.plot(magnetization, label='Magnetization')
plt.axhline(y=0, color='purple', ls='--')
plt.plot(energy, label='Energy')
plt.axhline(y=energy.mean(), color='red', ls='--')
plt.title(rf'Modelo de Ising: $c_H = {c_H}, \; \chi = {chi}$')
plt.legend(loc='upper right')
plt.show()