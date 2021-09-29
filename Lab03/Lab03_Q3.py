"""Created on 2021-09-29

@author: Gabriel Bailey

This code will 
"""
import numpy as np
import sympy as sp
from scipy.special import jv, kn
#from Lab02_IntegrationFunctions import J # Where the Bessel function is defined
import matplotlib.pyplot as plt

def f(x):
    return np.exp(-x**2)

def fd_derivative(f,x,h):
    return (f(x+h)-f(x))/(h)

def cd_derivative(f,x,h):
    return (f(x+h)-f(x-h))/(2*h)

x_val = 0.5
#Calculate the actual value 
x = sp.Symbol('x')
y = sp.exp(-x**2)
dx = sp.diff(y, x)
exact = dx.evalf(subs={x: x_val})

der_results = np.empty([3,17])
for i in range(17):
    h = 10**(-i)
    der_results[0][16-i] = h
    der_results[1][16-i] = fd_derivative(f,x_val,h)
    der_results[2][16-i] = abs(der_results[1][16-i] - exact)


print(der_results)
