# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np 
from numpy import empty,arange,exp,real,imag,pi
from numpy.fft import rfft,irfft

import matplotlib.pyplot as plt
#setting constants
L = 1 
d= 0.1
C = 1
sig = 0.3
N =  1000
x = np.arange(0,L,L/N)
v = 100
P = C*x*(L-x)/(L**2)
Q= (-(x-d)**2)/(2*(sig**2))
def dst(y):
    'function from dst package to perform 1D fourier sine-transform'
    N = len(y)
    y2 = empty(2*N,float)
    y2[0] = y2[N] = 0.0
    y2[1:N] = y[1:]
    y2[:N:-1] = -y[1:]
    a = -imag(rfft(y2))[:N]
    a[0] = 0.0

    return a


######################################################################
# 1D inverse DST Type-I

def idst(a):
    'function from dst package to perform 1D inverse-fourier sine-transform'
    N = len(a)
    c = empty(N+1,complex)
    c[0] = c[N] = 0.0
    c[1:N] = -1j*a[1:]
    y = irfft(c)[:N]
    y[0] = 0.0

    return y

phi_init = np.zeros(N)

psi_init = P*np.exp(Q)

phi_coeff = dst(phi_init)
psi_coeff = dst(psi_init)



def sol(t):
 'function that take in time value and calculates fourier coefficients for wave-function at that time'
 alpha = 0
 eta=0
 for k in range(1,N+1):
     omega= k*np.pi*t/L
     alpha += phi_coeff*np.cos(omega*t)
     eta += (psi_coeff*np.sin(omega*t))/omega

 total = alpha + eta
 return total 
sol_2 = idst(sol(0.002))
sol_4 = idst(sol(0.004))
sol_6 = idst(sol(0.006))
sol_12 = idst(sol(0.012))
sol_100 = idst(sol(0.1))

plt.plot(x,sol_2)
plt.title('Solution for vibration of string using Spectral method at 2ms')
plt.xlabel('distance along string(m)')
plt.ylabel('displacement of string(micro-meter)')
plt.show()

plt.plot(x,sol_4)
plt.title('Solution for vibration of string using Spectral method at 4ms')
plt.xlabel('distance along string(m)')
plt.ylabel('displacement of string(micro-meter)')
plt.show()

plt.plot(x,sol_6)
plt.title('Solution for vibration of string using Spectral method at 6ms')
plt.xlabel('distance along string(m)')
plt.ylabel('displacement of string(micro-meter)')
plt.show()

plt.plot(x,sol_12)
plt.title('Solution for vibration of string using Spectral method at 12ms')
plt.xlabel('distance along string(m)')
plt.ylabel('displacement of string(micro-meter)')
plt.show()

plt.plot(x,sol_100)
plt.title('Solution for vibration of string using Spectral method at 100ms')
plt.xlabel('distance along string(m)')
plt.ylabel('displacement of string(micro-meter)')
plt.show()
