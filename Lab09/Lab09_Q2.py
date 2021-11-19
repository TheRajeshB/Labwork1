# -*- coding: utf-8 -*-
'''Created on 2021-11-17

@author: Gabriel Bailey

This code will solve the time-dependant schrodinger equation with the Crank-Nicolson method, and will output a video file for one of the cases: b, c, or d.
'''
from numpy.random import rand
from numpy.linalg import solve
from numpy import dot,mean
from Lab09_functions import *
from time import time
import matplotlib.pyplot as plt
from scipy.constants import electron_mass, hbar
import numpy as np

from banded import banded

compute_psi = True
test = True

# Constants

L = 1e-8 #m
m = electron_mass #?
#print('m:',m)

sigma = L/25 #m
kappa = 500/L #m^-1

P = 1024
tau = 1e-18 #s
N = 3000

if test:
    P = 100
    tau = 1e-15
    N = 300
a = L/P #spacing of points
h = tau #rename
T = N*tau #s


# part a
x0 = L/5 #m

def wave(x):
    return np.exp(-(x-x0)**2/(4*sigma**2) + kappa*x*1j)


#psi0 = 1/symb_integrate(wave, -0.5*L, 0.5*L)
psi0 = 1/simp_integrate(wave, P*10, -0.5*L, 0.5*L)
print('psi0', psi0)



#Converts index to distance
def dis(i):
    return L*(i/P - 1/2)

#Converts distance to index
def ind(x):
    return int(np.round((x/L+1/2)*P))


# Initialize psi
psi = np.zeros((N,P),np.complex_)
for i in range(len(psi[0])):
    psi[0][i] = wave(dis(i))
psi[0] = psi0*psi[0]
# required functions

# Potential for a square well
def V(x): #Not used?
    if -L/2 < x < L/2:
        return 0
    else:
        return np.inf

# Discretized Hamiltonian
def HD(): # Not used, I just followed the textbook implementation
    return 0

# Computes the normalization diagnostic (integral of wavefunction?)
def compute_normalization(psi,i):
    return simp_integrate(lambda x : psi[i,int(x)], P-1, 0, P-1) * L/P
    #return simp_integrate(lambda x : psi[i,ind(x)], P-1, -0.5*L, 0.5*L)

#  Computes the energy
def compute_energy(psi):
    return 0

# Computes the expected position
def compute_expected_position(psi):
    return 0

# Calculate A & B matrices


a1 = 1 + h * hbar/(2*m*a**2) *1j
a2 = -h * hbar/(4*m*a**2) *1j
A = np.array([np.full(P,a2),np.full(P,a1),np.full(P,a2)])

b1 = 1 - h * hbar/(2*m*a**2) *1j
b2 = h * hbar/(4*m*a**2) *1j
B = np.diag(np.full(P,b1),0) + np.diag(np.full(P-1,b2),-1) + np.diag(np.full(P-1,b2),1)

def advance_timestep(psi,h,n):
    if n < N-1:
        v = B.dot(psi[n])
        psi[n+1] = banded(A,v,1,1)

if compute_psi:
    for i in range(N):
        advance_timestep(psi,h,i)

    np.savez('Q3_phi', psi=psi)
else:
    npzfile = np.load('Q3_phi.npz')
    psi = npzfile['psi']

#print(psi)
print('normalization 0',compute_normalization(psi,0))
print('normalization N',compute_normalization(psi,N-1))

#Plot the results as an animation
from matplotlib.animation import FuncAnimation, FFMpegWriter

fig, ax = plt.subplots()
ax.grid(which='both', axis='y')
ax.set_title('Probability') 
ax.set_xlabel('x')
ax.set_ylabel('Value')
xdata, ydata = [], []
line, = plt.plot([], [])
time_text = ax.text(0.1, 0.1, 'time = ', transform=ax.transAxes)

def init():
    ax.set_xlim(-L, L)
    ax.set_ylim(np.min(psi), np.max(psi))
    time_text.set_text('time = ')
    return line,

# Just plot static plots:
""" init()
x = linspace(-L/2, L/2, P)
y = psi[29,:]
line.set_data(x, y)
plt.plot(x,psi0*psi[0,:],label='bef')
plt.plot(x,psi0*psi[15,:])
plt.legend() """

def update(i):
    i = int(i)
    x = linspace(-L/2, L/2, P)
    y = psi[i,:]
    line.set_data(x, y)
    time_text.set_text('time = {:1.2e}'.format(i*h))
    return line, time_text

ani = FuncAnimation(fig, update, frames=N, interval=1,
                    init_func=init, blit=True)

# Save the animation (need ffmpeg)
""" f = 'Q3_vid.mp4' 
writervideo = FFMpegWriter(fps=60) 
ani.save(f, writer=writervideo) """

plt.show()