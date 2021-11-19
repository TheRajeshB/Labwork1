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
case = 'a'

# Constants

L = 1e-8 #m
m = electron_mass #?
#print('m:',m)

sigma = L/25 #m
kappa = 500/L #m^-1

P = 1024
tau = 1e-18 #s
N = 3000
V = 'function'

# Potential for a square well
def V_square(x):
    if -L/2 <= x <= L/2:
        return 0
    else:
        return np.inf

# Potential for a harmonic oscillator
omega = 3e15 # rad/s
def V_harm(x):
    if -L/2 <= x <= L/2:
        return 1/2*m*omega**2 * x**2
    else:
        return np.inf
    
# Potential for a Double well
V0 = 6e-17 # J
x1 = L/4
def V_double(x):
    if -L/2 <= x <= L/2:
        return V0*(x**2/x1**2-1)**2
    else:
        return np.inf

if case == 'a':
    V = V_square
    N = 3000
    x0 = L/5
elif case == 'c':
    V = V_harm
    N = 4000
    x0 = L/5
elif case == 'd':
    V = V_double
    N = 6000
    x0 = L/3

if test:
    P = 100
    tau = 1e-15
    N = 300
a = L/P #spacing of points
h = tau #rename
T = N*tau #s


# part a
x0 = L/5 #m

# The function for psi
def wave(x):
    return np.exp(-(x-x0)**2/(4*sigma**2) + kappa*x*1j)

# The function for psi* psi
def wave_abs(x):
    val = np.exp(-(x-x0)**2/(4*sigma**2) + kappa*x*1j)
    return val.conj()*val

#psi0 = 1/symb_integrate(wave, -0.5*L, 0.5*L)
psi0 = 1/np.sqrt(simp_integrate(wave_abs, P*10, -0.5*L, 0.5*L))

print('psi0', psi0)
print('psi0 res', psi0*simp_integrate(wave_abs, P*10, -0.5*L, 0.5*L))

#Converts index to distance
def dis(i):
    return L*(i/P - 1/2)

#Converts distance to index
def ind(x):
    return int((x/L+1/2)*P)
    return int(np.round((x/L+1/2)*P))


# required functions

# Discretized Hamiltonian
def HD(): # Not used, I just followed the textbook implementation
    return 0

# computes the probability of the particle being at a location given the value of the wave function there
def prob(psi):
    return psi.conj().T*psi

# Computes the normalization diagnostic (integral of wavefunction?)
def compute_normalization(psi,i):
    return simp_integrate(lambda x : prob(psi[i,int(x)]), P-1, 0, P-1)*L/P
    return simp_integrate(lambda x : psi[i,ind(x)].conj().T*psi[i,ind(x)], P, -0.5*L, 0.5*L-L/P)
    #return simp_integrate(lambda x : psi[i,ind(x)], P-1, -0.5*L, 0.5*L)


h1 = 1 - h * hbar/(2*m*a**2) *1j
h2 = h * hbar/(4*m*a**2) *1j
h1_arr = np.full(P,h1)
for i in range(len(h1_arr)):
    h1_arr[i] = h1_arr[i] + V(dis(i))
H = np.diag(h1_arr,0) + np.diag(np.full(P-1,h2),-1) + np.diag(np.full(P-1,h2),1)
# The internal part of the position integral
def ener(psi,x):
    return psi[i].conj().T* H *psi[i]
#  Computes the energy
def compute_energy(psi,i):
    return simp_integrate(lambda x : ener(psi,x), P-1, 0, P-1)*L/P

# The internal part of the position integral
def pos(psi,x):
    return psi.conj().T*x*psi
# Computes the expected position
def compute_expected_position(psi,i):
    return simp_integrate(lambda x : pos(psi[i,int(x)],x), P-1, 0, P-1)*L/P


# Initialize psi
psi = np.zeros((N,P),np.complex_)
for i in range(len(psi[0])):
    psi[0][i] = wave(dis(i))

print('norm           ',simp_integrate(wave_abs, P*10, -0.5*L, 0.5*L))
print('normalization-1',compute_normalization(psi,0))
psi[0] = psi0*psi[0]
print('normalization 0',compute_normalization(psi,0))


# Calculate A & B matrices


a1 = 1 + h * hbar/(2*m*a**2) *1j
a2 = -h * hbar/(4*m*a**2) *1j

b1 = 1 - h * hbar/(2*m*a**2) *1j
b2 = h * hbar/(4*m*a**2) *1j

a1_arr = np.full(P,a1)
b1_arr = np.full(P,b1)
for i in range(len(a1_arr)):
    a1_arr[i] = a1_arr[i] + 1j/(2*hbar)*V(dis(i))
    b1_arr[i] = b1_arr[i] - 1j/(2*hbar)*V(dis(i))

# We defined A with the diagonals as rows so that we can use banded.py
A = np.array([np.full(P,a2),a1_arr,np.full(P,a2)])
# Define B normally
B = np.diag(b1_arr,0) + np.diag(np.full(P-1,b2),-1) + np.diag(np.full(P-1,b2),1)
# print(A)
# print(B)

def advance_timestep(psi,h,n):
    if n < N-1:
        v = B.dot(psi[n])
        psi[n+1] = banded(A,v,1,1)

if compute_psi:
    for i in range(N):
        advance_timestep(psi,h,i)

    np.savez('Q3_phi_'+case, psi=psi)
else:
    npzfile = np.load('Q3_phi_'+case+'.npz')
    psi = npzfile['psi']

#print(psi)

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