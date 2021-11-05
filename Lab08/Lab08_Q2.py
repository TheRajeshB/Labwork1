# -*- coding: utf-8 -*-
"""
Created on Thu Nov  4 22:17:11 2021

@author: farzh
"""

from numpy import zeros, exp, linspace, array, copy
from pylab import plot, show, xlabel, ylabel, clf, pause,legend

# Constants
L = 1  # m
d = 0.1  # m
C = 1  # m/s
sigma = 0.3  # m
v = 100  # m/s
N = 1000  # number of spatial points
a = L / N  # spacing between points
h = 10 ** -6  # time step in seconds
k = h * v ** 2 / a ** 2
t_f = 2 * 10 ** -3


def psi_0(x):
    return C * (x * (L - x)) / L ** 2 * exp(-(x - d) ** 2 / (2 * sigma ** 2))


x_points = linspace(0, L, N + 1)
phi = zeros(N + 1, float)
psi = array(list(map(psi_0, x_points)), float)
results = []
t = 0.0
while t <= t_f:
    results.append(copy(phi))
    old_phi = copy(phi)
    old_psi = copy(psi)
    phi[1: N] = old_phi[1: N] + h * old_psi[1: N]
    psi[1: N] = old_psi[1: N] + k * (old_phi[2: N + 1] + old_phi[0: N - 1] - 2 * old_phi[1: N])
    t += h

print(len(results))
plot(x_points, results[0] * 10 ** 6, 'k',label='0')
plot(x_points, results[10] * 10 ** 6, 'b',label='10')
plot(x_points, results[100] * 10 ** 6, 'g',label='100')
plot(x_points, results[1100] * 10 ** 6)
legend()
xlabel("t (s)")
ylabel("$\phi$ (m)")
show()

