import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from matplotlib.widgets import Slider

def f(t, N, a, b, r):
    res = N * (r - a * (N - b)**2)

    return res

a = 1
b = 1
r = 2

N_0 = np.linspace(b - np.sqrt(r/a) - 1, b + np.sqrt(r/a) + 1, 10)
t_span = [0, 3]
t_eval = np.linspace(*t_span, 100)

result = solve_ivp(f, t_span, N_0, t_eval=t_eval, args=[a, b, r])

t = result.t
N = result.y.T

fig,ax = plt.subplots()
fig.subplots_adjust(bottom=0.3)

lines = ax.plot(t, N)
ax.set_xlabel("t")
ax.set_ylabel("N(t)")
ax.grid()

# Sliders
axa = fig.add_axes((0.25, 0.1, 0.65, 0.03))
axb = fig.add_axes((0.25, 0.05, 0.65, 0.03))
axr = fig.add_axes((0.25, 0.0, 0.65, 0.03))

a_slider = Slider(ax=axa, label="a", valmin=1, valmax=10, valinit=1)
b_slider = Slider(ax=axb, label="b", valmin=0, valmax=10, valinit=1)
r_slider = Slider(ax=axr, label="r", valmin=0, valmax=10, valinit=1)

def update(val):

    # Obtener nuevos parámetros
    a = a_slider.val
    b = b_slider.val
    r = r_slider.val

    # Actualizar condiciones iniciales
    N_0 = np.linspace(
        b - np.sqrt(r / a) - 1,
        b + np.sqrt(r / a) + 1,
        10
    )

    # Resolver de nuevo
    result = solve_ivp(
        f,
        t_span,
        N_0,
        t_eval=t_eval,
        args=(a, b, r)
    )

    t = result.t
    N = result.y.T

    # Actualizar cada curva
    for i, line in enumerate(lines):
        line.set_data(t, N[:, i])

    # Actualizar límites del gráfico
    ax.relim()
    ax.autoscale_view()

    fig.canvas.draw_idle()

a_slider.on_changed(update)
b_slider.on_changed(update)
r_slider.on_changed(update)

plt.show()