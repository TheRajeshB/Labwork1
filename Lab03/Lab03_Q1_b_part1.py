# -*- coding: utf-8 -*-
"""
Created on Tue Sep 28 23:17:56 2021

@author: farzh
"""

from functions_Lab03 import simp_integrate, trap_integrate , gaussxwab
from math import *
import matplotlib.pyplot as plt
import numpy as np
import scipy.special as sc
import matplotlib.pyplot as plt

def u(x,lamb,z):
    ''' function to calculate values of u
        x = x-coordinate of position
        z = z coordinate of postion
        lamb = wavelength'''
    return x*((2/(lamb*z))**0.5)





def f(t):
    'function for integrand of C(u)'
    return cos(0.5*pi*(t**2))

def g(t):
    'function of integrand of S(u)'
    return sin(0.5*pi*(t**2))

def C(u,N,a,x,lamb,z):
    '''function to calculate C(u)
       a = start point of integration
       u = function defining values of b(end point of integration)
       N = number of slices
       u = function defined earlier
       x = x-coordinate of position
       z = z coordinate of postion
       lamb = wavelength'''
    C_values = []
    for i in u(x,lamb,z):
     b = i   
     x1,w1 = gaussxwab(N,a,b)
     s1 = 0.0
     for k in range(N):
        s1 += w1[k]*f(x1[k])
     C_values.append(s1)
    return C_values


def S(u,N,a,x,lamb,z):
    '''function to calculate S(u)
       a = start point of integration
       u = function defining values of b(end point of integration)
       N = number of slices
       u = function defined earlier
       x = x-coordinate of position
       z = z coordinate of postion
       lamb = wavelength'''
       
    S_values= []
    for i in u(x,lamb,z):
           b = i
           x2,w2 = gaussxwab(N,a,b)
           s2 = 0.0
           for k in range(N):
               s2 += w2[k]*g(x2[k])
           S_values.append(s2)
    return S_values       

p =np.array([ u(np.linspace(-5,5),1,3)]) #        
q = np.array([C(u,50,0.0,np.linspace(-5,5),1,3)]) #assigning name for results of C(u), using Gaussian, for simplicity later
r =np.array ([S(u,50,0.0,np.linspace(-5,5),1,3)])#assigning name for results of S(u), using Gaussian, for simplicity later

S,C = sc.fresnel(p) #reslts for S(u) and C(u) from scipy
I_ratio_Gauss = (((1+2*q)**2) + ((1+2*r)**2))/8 #Intensity ratio using Gaussian
I_ratio_scipy = (((1+2*C)**2) + ((1+2*S)**2))/8 #Intensity ratio using scipy
I_err = abs((I_ratio_scipy - I_ratio_Gauss)/I_ratio_scipy) #relative error of results from Gaussian with respect to results from Scipy


plt.errorbar(p,C,marker='o',markersize=3, label = 'result for C(u) from scipy')
plt.errorbar(p,q,label='result from C(u) from Gaussian')
plt.xlabel('values of u')
plt.ylabel('values for C(u)' )
plt.title('Comparison of integration results using Scipy and Gaussian for C(u)')
plt.legend()
plt.show()


plt.errorbar(p,S,marker='o',markersize=3, label = 'result for S(u) from scipy')
plt.errorbar(p,r,label='result from S(u) from Gaussian')
plt.xlabel('values of u')
plt.ylabel('values for S(u)' )
plt.title('Comparison of integration results using Scipy and Gaussian for S(u)')
plt.legend()
plt.show()

plt.errorbar(p,(((1+2*q)**2) + ((1+2*r)**2))/8,marker='o',label='Reslts for I/I0 using Gaussian')
plt.errorbar(p,(((1+2*C)**2) + ((1+2*S)**2))/8, label='results for I/I0 using scipy')
plt.xlabel('values of u')
plt.ylabel('values for I/I0' )
plt.title('Comparison of results for I/I0 using Gaussian and Scipy')
plt.legend()
plt.show()

plt.errorbar(p,I_err,marker='o',label='relative errr')
plt.xlabel('values of u')
plt.ylabel('relative error for values for I/I0' )
plt.title('Relative error of results from Gaussian with respect to results from Scipy')
plt.legend()
plt.show()

