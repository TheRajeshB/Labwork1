# SolveLinear.py
# Python module for PHY407
# Paul Kushner, 2015-09-26
# Modifications by Nicolas Grisouard, 2018-09-26
# This module contains useful routines for solving linear systems of equations.
# Based on gausselim.py from Newman

import numpy as np
from scipy.constants import pi,c
from numpy import ones,empty,copy,cos,tan,linspace

# From Lab 3:

def gaussxw(N):

    # Initial approximation to roots of the Legendre polynomial
    a = linspace(3,4*N-1,N)/(4*N+2)
    x = cos(pi*a+1/(8*N*N*tan(a)))

    # Find roots using Newton's method
    epsilon = 1e-15
    delta = 1.0
    while delta>epsilon:
        p0 = ones(N,float)
        p1 = copy(x)
        for k in range(1,N):
            p0,p1 = p1,((2*k+1)*x*p1-k*p0)/(k+1)
        dp = (N+1)*(p0-x*p1)/(1-x*x)
        dx = p1/dp
        x -= dx
        delta = max(abs(dx))

    # Calculate the weights
    w = 2*(N+1)*(N+1)/(N*N*(1-x*x)*dp*dp)

    return x,w

def gaussxwab(N,a,b):
    x,w = gaussxw(N)
    return 0.5*(b-a)*x+0.5*(b+a), 0.5*(b-a)*w

def gaussxwab_convert(x,w,a,b):
    return 0.5*(b-a)*x+0.5*(b+a), 0.5*(b-a)*w

# A function to calculate the velocity of a spring at position x with initial (0 velocity) position x0
def g(x,x0,k,m):
    return c*( k*(x0**2-x**2)*(2*m*c**2 + k*(x0**2-x**2)/2) / (2*(m*c**2 + k*(x0**2-x**2)/2)**2) )**(1/2)

# Calculates T using Gaussian quadrature for the integral
def T_details(x0, k, m, N):
    x,w = gaussxwab(N,0,x0)
    gk, wgk = empty(N), empty(N)
    sum = 0.0
    for i in range(N):
        gk[i] = 4/g(x[i],x0,k,m)
        wgk[i] = w[i]*gk[i]
        sum += wgk[i]
    return sum, x, gk, wgk




