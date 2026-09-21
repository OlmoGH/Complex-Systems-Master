import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import RegularPolyCollection
from matplotlib.animation import FuncAnimation
from matplotlib.widgets import Slider

def montecarlo_step(grid, T, J=-1):

    L = grid.shape[0]

    i, j = np.random.randint(0, L, 2)

    spin = grid[i, j]

    if i % 2 == 0:
        neighbors = [
            (i, j-1),
            (i, j+1),
            (i-1, j-1),
            (i-1, j),
            (i+1, j-1),
            (i+1, j),
        ]
    else:
        neighbors = [
            (i, j-1),
            (i, j+1),
            (i-1, j),
            (i-1, j+1),
            (i+1, j),
            (i+1, j+1),
        ]

    neighbor_sum = sum(
        grid[ii % L, jj % L]
        for ii, jj in neighbors
    )

    dE = 2 * J * spin * neighbor_sum

    if np.random.rand() < np.exp(dE / T):
        grid[i, j] *= -1

L = 100
T = 5
grid = np.random.choice([-1, 1], size=(L, L))

# Creation of hexagonal points
x = [np.sqrt(3) * (j + 0.5 * (i % 2)) for i in range(L) for j in range(L)]
y = [1.5 * i for i in range(L) for j in range(L)]



fig, ax = plt.subplots(figsize=(7, 7))
collection = RegularPolyCollection(numsides=6, 
                                   rotation=0, 
                                   sizes=(10,), 
                                   offsets=np.column_stack([x, y]), 
                                   cmap='binary',
                                   offset_transform=ax.transData,
                                   edgecolor='None')
collection.set_array(grid.ravel())
collection.set_clim(-1, 1)
ax.add_collection(collection)
ax.set_aspect("equal")
ax.autoscale()
ax.axis("off")
fig.subplots_adjust(bottom=0.25)
T_ax = fig.add_axes([0.25, 0.1, 0.65, 0.03])
T_slider = Slider(T_ax, 'T', valmin=0.1, valmax=5, valinit=1)


def update(frame):
    montecarlo_step(grid, T_slider.val)

    collection.set_array(grid.ravel())

    return collection, 



animation = FuncAnimation(fig, func=update, frames=None, cache_frame_data=False, interval=20)

plt.show()