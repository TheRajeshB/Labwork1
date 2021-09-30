# -*- coding: utf-8 -*-
"""
Created on Wed Sep 29 23:46:57 2021

@author: farzh
"""
from functions_Lab03 import simp_integrate, trap_integrate , gaussxwab
from math import *
import matplotlib.pyplot as plt
import numpy as np
import scipy.special as sc

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
    '''function to calculate C(u) and store in list
       a = start point of integration
       u = function defining values of b(end point of integration)
       N = number of slices
       u = function defined earlier
       x = x-coordinate of position
       z = z coordinate of postion
       lamb = wavelength'''
    C_values = []
    for i in u(x,lamb,z):
     b = i   #end points of integration
     x1,w1 = gaussxwab(N,a,b)
     s1 = 0.0
     for k in range(N): #loop for Gaussian integration
        s1 += w1[k]*f(x1[k])
     C_values.append(s1)
    return C_values

def S(u,N,a,x,lamb,z):
    '''function to calculate S(u) and store in list
       a = start point of integration
       u = function defining values of b(end point of integration)
       N = number of slices
       u = function defined earlier
       x = x-coordinate of position
       z = z coordinate of postion
       lamb = wavelength'''
       
    S_values= []
    for i in u(x,lamb,z):
           b = i #end points of integration
           x2,w2 = gaussxwab(N,a,b)
           s2 = 0.0
           for k in range(N): #loop for Gaussian integration
               s2 += w2[k]*g(x2[k])
           S_values.append(s2)
    return S_values 

a = np.array([C(u,50,0.0,np.linspace(-3,10,10),1,np.linspace(1,5,10))])
b = np.array([S(u,50,0.0,np.linspace(-3,10,10),1,np.linspace(1,5,10))])
I1 = (((2*a + 1)**2) + ((2*b + 1)**2))/8
list = []
for i in I1 :
    for j in i:
        list.append(j)
        
print(list)
x = np.linspace(-3,10,5)
z = np.linspace(1,5,5)
z,x = np.meshgrid(z,x)
