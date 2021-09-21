"""Created on 2021-09-21

@author: Gabriel Bailey

This code will calculate definite integrals using Trapezoidal and Simpsons's rules for integration.
"""

import numpy as np
import matplotlib.pyplot as plt
from time import time

import sympy # the symbolic math package

def symb_integrate(f, a, b): # symbolic integration to test errors
    #Directly from the lecture slides:
    xs = sympy.Symbol('xs', real=True) # the variable of integration
    return(sympy.integrate(f(xs), (xs, a, b)))

def trap_integrate(f, N, a, b):
    #Basically directly from the lecture slides:
    h = (b-a)/N
    area = 0.5*f(a) + 0.5*f(b) # the end bits
    for k in range(1,N): # adding the interior bits
        area += f(a+k*h)
    return h*area

def simp_integrate(f, N, a, b):
    #Based on the lecture slides:
    h = (b-a)/N
    area = 0
    if N // 2 != 0: #if N is odd calc the last slice with trapezoidal rule
        area += f(a+(N-1)*h)
        N -= 1
    
    A_odd = 0
    for k in range(1,N,2): # for the odd terms
        #print(A_odd*h,h,k, N)
        A_odd += f(a+k*h)

    A_even = 0
    for k in range(2,N,2): # for the odd terms
        #print(A_even*h,h,k, N)
        A_even += f(a+k*h)

    area += h/3 * (f(a) + f(b) + 4*A_odd + 2*A_even)
    return area

def f(x):
    return(4/(1+x**2))

#ii

print("Real value (via sympy):    ",float(symb_integrate(f,0,1)))
print("Trapezoidal integration:   ",trap_integrate(f,10,0,1))
print("Simpson's rule integration:",simp_integrate(f,10,0,1))

#iii


val = symb_integrate(f,0,1)


n = 10
trap_val = trap_integrate(f,n,0,1)
while float(abs(trap_val-val)) > 10**-9: #Multiplicative loop to get to general vicinity fast
    n += int(n*0.1)
    trap_val = trap_integrate(f,n,0,1)

n = int(n /1.1) #backtrack one
trap_val = trap_integrate(f,n,0,1)
while float(abs(trap_val-val)) > 10**-9: #One by one to get "exact" point
    #print(float(abs(trap_val-val)))
    n += 1
    trap_val = trap_integrate(f,n,0,1)

print("Trap N required to reduce error below O(10^-9):", n)
print("Trap error:", float(abs(trap_val-val)))
trap_start = time()
for i in range(1000):
    trap_integrate(f,n,0,1)
trap_end = time()
print("Trap Time Taken:",(trap_end-trap_start)/1000, "s")
'''
n = 4
simp_val = simp_integrate(f,n,0,1)
while abs(simp_val-val) > 10**-9:
    print(abs(simp_val-val))
    n += 1
    simp_val = simp_integrate(f,n,0,1)

print("Simp N requied to reduce error below O(10^-9):", n)
print("Simp error:", abs(simp_val-val))
'''
'''
# Plot a time comparision of our histogram function to numpy's histogram function with a logarithmic y axis.
plt.plot(Ns, times, color = 'b', label="My histogram")
plt.plot(Ns, nptimes, color = 'r', label = "numpy's histogram")
plt.legend()
plt.xlabel('N')
plt.ylabel('Time taken (s)')
plt.title('Histogram Function Time Comparision')
plt.yscale("log")
#plt.xscale("log")
plt.show()
'''