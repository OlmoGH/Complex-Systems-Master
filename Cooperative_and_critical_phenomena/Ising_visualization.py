import numpy as np
import matplotlib.pyplot as plt
import os
import scipy.signal as sc

directorio = os.path.dirname(os.path.abspath(__file__))

data = np.loadtxt(directorio + "/out/Data.txt", skiprows=1, delimiter='\t').T
magnetization = np.loadtxt(directorio + "/out/Correlation.txt")

def calcular_tiempo_correlacion(magnetization, umbral_corte=0.05):
    M = np.array(magnetization)
    N = len(M)
    
    # 1. Aislar las fluctuaciones restando la media térmica
    M_fluct = M - np.mean(M)
    
    # 2. Calcular la autocorrelación directa
    # mode='full' devuelve 2N-1 elementos, seleccionamos la mitad derecha (tau >= 0)
    corr = np.correlate(M_fluct, M_fluct, mode='full')[N-1:]
    
    # 3. Corregir el sesgo estadístico por la pérdida de pares solapados
    lags = np.arange(N)
    corr_insesgada = corr / (N - lags)
    
    # 4. Normalizar para que C(0) sea exactamente 1.0
    C_tau = corr_insesgada / corr_insesgada[0]
    
    # 5. Truncar la suma para evitar que el ruido estadístico de la cola domine
    # Se corta cuando la correlación cae por debajo del umbral o cruza el cero
    cruces = np.where(C_tau < umbral_corte)[0]
    tau_max = cruces[0] if len(cruces) > 0 else N
    
    # 6. Calcular el tiempo de correlación integrado
    tau_int = 0.5 + np.sum(C_tau[1:tau_max])
    
    return tau_int, C_tau[:tau_max]

# Ejemplo de uso:
# magnetization = [...] # Tu serie temporal extraída de la simulación
tau, C = calcular_tiempo_correlacion(magnetization)
print(f"Tiempo de correlación integrado: {tau:.2f} pasos")

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