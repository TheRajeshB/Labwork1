"""Created on 2021-09-29

@author: Gabriel Bailey

This code will calculate derivatives with the forward difference and centered difference methods as well as their errors over a variety of values for h. It will output a plot of the errors and print the results for forward difference in a csv table.
"""
import numpy as np
#import sympy as sp
from scipy.special import jv, kn
#from Lab02_IntegrationFunctions import J # Where the Bessel function is defined
import matplotlib.pyplot as plt

#The function in question
def f(x):
    return np.exp(-x**2)

# Analytic derivative of f
def dfdx(x):
    return -2*x*np.exp(-x**2)

# Forward difference method
def fd_derivative(f,x,h):
    return (f(x+h)-f(x))/(h)

# Central difference method
def cd_derivative(f,x,h):
    return (f(x+h)-f(x-h))/(2*h)

# Set the x value we're interested in and calculate the exact derivative (with roundoff errors obv)
x_val = 0.5
exact = dfdx(x_val)

#Calculate the actual value generally using sympy
'''x = sp.Symbol('x')
y = sp.exp(-x**2)
dx = sp.diff(y, x)
exact = dx.evalf(subs={x: x_val})'''

# Initialize data storage lists for h and errors
hs = np.empty(17)
fd_diff = np.empty(17)
cd_diff = np.empty(17)

print("h,fd value,fd error")
for i in range(0,17):
    # Calculate the derivatives and errors with forward difference and centered difference, then print 3b table
    h = 10**(i-16)
    hs[i] = h
    fd_diff[i] = abs(fd_derivative(f,x_val,h) - exact)
    cd_diff[i] = abs(cd_derivative(f,x_val,h) - exact)
    print(hs[i],',',fd_derivative(f,x_val,h),',',fd_diff[i])


#Plot errors of the two methods on the same plot
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
