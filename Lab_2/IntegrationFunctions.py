"""Created on 2021-09-21

@author: Gabriel Bailey

This code will define functions that calculate definite integrals using Symbolic maths, Trapezoidal rule and Simpsons's rules for integration.
"""
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