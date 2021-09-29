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

hs = np.empty(17)
fd_val = np.empty(17)
fd_diff = np.empty(17)
cd_diff = np.empty(17)

print("h,fd value,fd error")
for i in range(0,17): #i = 16-i
    h = 10**(i-16)
    hs[i] = h
    fd_diff[i] = abs(fd_derivative(f,x_val,h) - exact)
    cd_diff[i] = abs(cd_derivative(f,x_val,h) - exact)
    print(hs[i],',',fd_derivative(f,x_val,h),',',fd_diff[i])


#Plot absolute difference from scipy Bessel function values
plt.figure()
plt.plot(hs, fd_diff, label = "Forward Difference Error")
plt.plot(hs, cd_diff, label = "Central Difference Error")
plt.xscale('log')
plt.yscale('log')
plt.legend()
plt.xlabel('h')
plt.ylabel('Absolute Error')
plt.title('Errors for Derivative of e^(-x^2) at x=0.5')  

plt.show()
