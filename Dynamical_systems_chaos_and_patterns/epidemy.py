import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import root

def f(a, b, u):

    return a - b * u - np.exp(-u) 

L = 10
b_range = np.linspace(0, 10, L)
a_range = np.linspace(1, 10, L)

X, Y = np.meshgrid(a_range, b_range, indexing='ij')

M = np.stack((X, Y), axis=-1)

grid = np.zeros((L, L))
