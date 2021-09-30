# -*- coding: utf-8 -*-
"""
Created on Wed Sep 29 20:50:19 2021

@author: farzh
"""
from functions_Lab03 import simp_integrate, trap_integrate , gaussxwab
from math import *
import matplotlib.pyplot as plt
import numpy as np
import scipy.special as sc

#from previous question we notice signal becomes smooth around u >0
#the corresponding x values show that this is for x>0


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
N_list=[]#list to store values of N
max_delta_list =[] #list to store values of max of relative difference for each N value
for j in range(3,51):
    
    N=j #values for N
    N_list.append(j)
    S_Gauss = np.array([S(u,N,0.0,np.linspace(0,5),1,3)]) #values for S(u) using Gaussian integration
    C_Gauss = np.array([C(u,N,0.0,np.linspace(0,5),1,3)]) #values for C(u) using Gaussian integration
    S_scipy,C_scipy = sc.fresnel(u(np.linspace(0,5),1,3)) #values for S(u) and C(u) using Scipy
    I_ratio_Gauss = (((1+2*C_Gauss)**2) + ((1+2*S_Gauss)**2))/8 #Intensity ratio using Gaussian
    I_ratio_scipy = (((1+2*C_scipy)**2) + ((1+2*S_scipy)**2))/8 #Intensity ratio using scipy
    
    
    relative_diff = abs((I_ratio_scipy-I_ratio_Gauss)/I_ratio_scipy) #relative error of results from Gaussian with respect to results from Scipy
    
    max_delta_list.append(np.max(relative_diff)) #add max of relative difference to list
    

plt.errorbar(N_list,max_delta_list,marker='o',label='maximum of relative diff for corresponding N value')
plt.xlabel('value of N')
plt.ylabel('maximum of relative difference')
plt.title('Graph showing maximum of relative error for range of N-values ')
plt.legend()
plt.show()
