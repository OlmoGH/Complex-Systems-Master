import numpy as np
import matplotlib.pyplot as plt
import os
import scipy.signal as sc

directorio = os.path.dirname(os.path.abspath(__file__))

data = np.loadtxt(directorio + "/out/Data.txt", skiprows=1, delimiter='\t').T

T = data[0]
E = data[1]
M = data[2]
C_H = data[3]
chi = data[4]

fig, ax = plt.subplots(ncols=2, nrows=2)

ax[0, 0].plot(T, E, linestyle='', marker='.')
ax[0, 0].set(title="Energy")

ax[0, 1].plot(T, M, linestyle='', marker='.')
ax[0, 1].set(title="Magnetization")

ax[1, 0].plot(T, C_H, linestyle='', marker='.')
ax[1, 0].set(title=r"$C_H$")

ax[1, 1].plot(T, chi, linestyle='', marker='.')
ax[1, 1].set(title=r"$\chi$")

plt.tight_layout()
plt.show()