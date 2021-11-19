# -*- coding: utf-8 -*-
'''Created on 2021-11-17

@author: Gabriel Bailey

This code will solve the time-dependant schrodinger equation with the Crank-Nicolson method, and will output a video file for one of the cases: b, c, or d. It will also save data if needed.
'''
from math import inf
from Lab09_functions import *
import matplotlib.pyplot as plt
from scipy.constants import electron_mass, hbar
import numpy as np

from banded import banded

# Options
compute_psi = True
compute_X = False
compute_en = False
test = False # For quick testing
case = 'd'
save_vid = False

# Constants for this problem
L = 1e-8 #m
m = electron_mass #kg

sigma = L/25 #m
kappa = 500/L #m^-1

P = 1024
tau = 1e-18 #s
N = 3000
V = 'function'
x0 = L/5 #m

# Potential functions:

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

# Differentiate between cases
if case == 'b':
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

# required functions

# Discretized Hamiltonian
A1 = - hbar**2/(2*m*a**2)
Bp = np.zeros(P,np.complex_)
for i in range(P):
    Bp[i] = V(dis(i)) + -2*A1
HD = np.diag(Bp,0) + np.diag(np.full(P-1,A1),-1) + np.diag(np.full(P-1,A1),1)

# computes the probability of the particle being at a location given the value of the wave function there
def prob(psi):
    return psi.conj()*psi

# Computes the normalization diagnostic (integral of wavefunction?)
def compute_normalization(psi,i):
    return simp_integrate(lambda x : prob(psi[i,int(x)]), P-1, 0, P-1)*L/P
    #return simp_integrate(lambda x : psi[i,ind(x)].conj().T*psi[i,ind(x)], P, -0.5*L, 0.5*L-L/P)
    #return simp_integrate(lambda x : psi[i,ind(x)], P-1, -0.5*L, 0.5*L)

#  Computes the energy
def ener(psi,i):

    #print('once',psi.conj().T.shape, HD.shape, psi.shape)
    print('What is going on with these matrix shapes, v_T x H = v x A')
    print('??? starter:',psi.conj().T.shape, psi.conj().shape)
    print('??? result :',(psi.conj().T*HD).shape, (psi.conj()*HD).shape)
    print('final',(psi.conj().T* HD* psi).shape)
    return (psi.conj().T* HD *psi)[i]
    
def compute_energy(psi,i):
    en = ener(psi[i][np.newaxis],i)
    return simp_integrate(lambda x: en[int(x)], P-1, 0, P-1)*L/P

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
# Normalize psi
psi[0] = psi0*psi[0]

# Print some info
print('normalization 0',compute_normalization(psi,0))
print('position 0',compute_expected_position(psi, 0))
#print('energy 0',compute_energy(psi,0))

# Calculate A & B matrices from textbook

a1 = 1 + h * hbar/(2*m*a**2) *1j
a2 = -h * hbar/(4*m*a**2) *1j

b1 = 1 - h * hbar/(2*m*a**2) *1j
b2 = h * hbar/(4*m*a**2) *1j

a1_arr = np.full(P,a1)
b1_arr = np.full(P,b1)
for i in range(len(a1_arr)):
    a1_arr[i] = a1_arr[i] + h*1j/(2*hbar)*V(dis(i))
    b1_arr[i] = b1_arr[i] - h*1j/(2*hbar)*V(dis(i))

# We defined A with the diagonals as rows so that we can use banded.py
A = np.array([np.full(P,a2),a1_arr,np.full(P,a2)])
# Define B normally
B = np.diag(b1_arr,0) + np.diag(np.full(P-1,b2),-1) + np.diag(np.full(P-1,b2),1)

def advance_timestep(psi,h,n):
    if n < N-1:
        v = B.dot(psi[n])
        psi[n+1] = banded(A,v,1,1)

if compute_psi:
    for i in range(N):
        advance_timestep(psi,h,i)

    np.savez('Q2_psi_'+case, psi=psi)
else:
    npzfile = np.load('Q2_psi_'+case+'.npz')
    psi = npzfile['psi']

X = np.zeros(N,np.complex_)
# Uncomment to try to calculate X
""" if compute_X:
    for i in range(N):
        X[i] = compute_expected_position(psi,i)

    np.savez('Q2_X_'+case, X=X)
else:
    npzfile = np.load('Q2_X_'+case+'.npz')
    X = npzfile['X'] """

energy = np.zeros(N,np.complex_)
# Uncomment to try to calculate energy
""" if compute_en:
    for i in range(N):
        X[i] = compute_energy(psi,i)

    np.savez('Q2_en_'+case, en=en)
else:
    npzfile = np.load('Q2_en_'+case+'.npz')
    energy = npzfile['en'] """
#print(psi)

print('normalization N',compute_normalization(psi,N-1))
#print('energy N',compute_energy(psi,N-1))
print('position N',X[N-1])
#Plot the results as an animation
from matplotlib.animation import FuncAnimation, FFMpegWriter

# plot normalization over time
""" norm = np.zeros(N)
for i in range(N):
        norm[i] = compute_normalization(psi,i)
figdiag, axdiag = plt.subplots()
axdiag.grid(which='both', axis='y')
axdiag.set_title('Part '+case+' Normalization') 
axdiag.set_xlabel('Time (s)')
axdiag.set_ylabel('Normalization')
line, = plt.plot(np.linspace(0,T,N), norm)
plt.show() """

# Create Animation
fig, ax = plt.subplots()
ax.grid(which='both', axis='y')
ax.set_title('Part '+case+' Phi (Real Part)') 
ax.set_xlabel('Position (m)')
ax.set_ylabel('Real Value')
xdata, ydata = [], []
line, = plt.plot([], [])
time_text = ax.text(0.1, 0.1, 'time = ', transform=ax.transAxes)
#vline = ax.axvline(X[0])

def init():
    ax.set_xlim(-6/10*L, 6/10*L)
    ax.set_ylim(np.min(psi), np.max(psi))
    time_text.set_text('time = ')
    return line,

def update(i):
    i = int(i)
    x = linspace(-L/2, L/2, P)
    y = psi[i,:]
    line.set_data(x, y)
    time_text.set_text('time = {:1.2e}'.format(i*h))
    #vline.set_data( [i, i], [-inf, inf])
    return line, time_text

ani = FuncAnimation(fig, update, frames=N, interval=1,
                    init_func=init, blit=True)

# Save the animation (need ffmpeg)
if save_vid:
    f = 'Q2_vid_'+case+'.mp4' 
    writervideo = FFMpegWriter(fps=60) 
    ani.save(f, writer=writervideo)

plt.show()