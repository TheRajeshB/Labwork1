"""Created on 2021-09-29

@author: Gabriel Bailey

This code will 
"""

import numpy as np
#import sympy as sp
from scipy.constants import c
from functions_Lab03 import gaussxwab
import matplotlib.pyplot as plt

def g(x,x0,k,m):
    return c*( k*(x0**2-x**2)*(2*m*c**2+k*(x0**2-x**2)/2) / (2*(m*c**2+k*(x0**2-x**2)/2)**2) )**(1/2)

def T(x0,k,m, N):
    x,w = gaussxwab(N,0,x0)
    s = 0.0
    for i in range(N):
        s += w[i]*1/g(x[i],x0,k,m)
    return 4*s

x0 = 0.01   # m
m = 1       # kg
k = 12      # N/m
T8 = T(x0,k,m,8)
T16 = T(x0,k,m,8)
print(T8,T16)
err = abs((T16-T8)/T8)
print(err)
