# SolveLinear.py
# Python module for PHY407
# Paul Kushner, 2015-09-26
# Modifications by Nicolas Grisouard, 2018-09-26
# This module contains useful routines for solving linear systems of equations.
# Based on gausselim.py from Newman
from numpy import empty
# The following will be useful for partial pivoting
from numpy import empty, copy
from numpy import array

#FUNCTIONS FOR LAB 4

def GaussElim(A_in, v_in):
    """Implement Gaussian Elimination. This should be non-destructive for input
    arrays, so we will copy A and v to
    temporary variables
    IN:
    A_in, the matrix to pivot and triangularize
    v_in, the RHS vector
    OUT:
    x, the vector solution of A_in x = v_in """
    # copy A and v to temporary variables using copy command
    A = copy(A_in)
    v = copy(v_in)
    N = len(v)
    

    
    for m in range(N):
        # Divide by the diagonal element
        div = A[m, m]
        A[m, :] /= div
        v[m] /= div

        # Now subtract from the lower rows
        for i in range(m+1, N):
            mult = A[i, m]
            A[i, :] -= mult*A[m, :]
            v[i] -= mult*v[m]

    # Backsubstitution
    # create an array of the same type as the input array
    x = empty(N, dtype=v.dtype)
    for m in range(N-1, -1, -1):
        x[m] = v[m]
        for i in range(m+1, N):
            x[m] -= A[m, i]*x[i]
    return x


def PartialPivot(A_in, v_in):
    """ In this function, code the partial pivot (see Newman p. 222) """
    'hello'
    A = copy(A_in)
    v = copy(v_in)
    N = len(v)
    #performing partial pivoting 
    
    for i in range(N): 
        for j in range(i+1,N):
         if abs(A[i,i]) < abs(A[j,i]): #comparing mangitudes of each diagonal element to magnitudes of all elements below them
            A[i, :], A[j, :] = copy(A[j, :]), copy(A[i, :]) #swapping rows so that each diagonal element has highest magnitude in its column
            v[i], v[j] = copy(v[j]), copy(v[i]) #performing same swap for corresponding rows in v
    
    #now matrix and vectors have been partial pivoted, perform Gaussian elimination:
    for m in range(N):
        # Divide by the diagonal element
        div = A[m, m]
        A[m, :] /= div
        v[m] /= div

        # Now subtract from the lower rows
        for i in range(m+1, N):
            mult = A[i, m]
            A[i, :] -= mult*A[m, :]
            v[i] -= mult*v[m]

    # Backsubstitution
    # create an array of the same type as the input array
    x = empty(N, dtype=v.dtype)
    for m in range(N-1, -1, -1):
        x[m] = v[m]
        for i in range(m+1, N):
            x[m] -= A[m, i]*x[i]
    return x

#FUNCTIONS FROM LAB 3

"""
This code will define functions that calculate definite integrals using Symbolic maths, Trapezoidal rule and Simpsons's rules for integration.
"""
import sympy # the symbolic math package
import numpy as np
from scipy.constants import pi
from numpy import ones,copy,cos,tan,pi,linspace

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
    if N % 2 == 1: #if N is odd calc the last slice with trapezoidal rule
        area += h*f(a+(N-1)*h)
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

#Bessel functions:
def f_inner(phi,x,n):
    return(np.cos(n*phi-x*np.sin(phi)))

# nth order Bessel function
def J(x,n):
    return(1/pi*simp_integrate(lambda phi : f_inner(phi,x,n),1000,0,pi))

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
            





