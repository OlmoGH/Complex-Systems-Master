import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("Simulation_methods/Errores.txt", delimiter='\t', skiprows=1)

M = data[:, 0]
err = data[:, 1]

x = np.linspace(10, 1e9, 1000)
a, b = np.polyfit(np.log(M), np.log(err), 1)

plt.loglog(x, np.exp(b) * x ** a, 'r-', label="Linear fit")
plt.loglog(x, np.sqrt(np.pi/4 * (1- np.pi / 4) / x), 'g--', label="Theoretical error")
plt.loglog(M, err, 'b.', label="Data")
plt.title(rf"$\alpha$ = {a: .3f}")
plt.legend()
plt.show()