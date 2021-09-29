# -*- coding: utf-8 -*-
"""
Created on Wed Sep 29 13:51:19 2021

@author: farzh
"""

from functions_Lab03 import simp_integrate, trap_integrate , gaussxwab
from math import *
import matplotlib.pyplot as plt


#defining function to be integrated
def f(x):
    return 4/(1+x**2)

simp_values = []
trap_values = []
Gauss_values = []
Gauss_values_2N=[]
N_values = []
a=0
b=1

for i in range(3,12,1): #N is even numbers between from 2 to 2048 due to Simpson integration 
    
    N = 2**i
    N_values.append(N)
    
    simp_values.append(simp_integrate(f,N,a,b)) 
    trap_values.append(trap_integrate(f,N,0,1))
    x,w = gaussxwab(N,a,b)
    s = 0.0
    s2= 0.0
    for k in range(N):
        s += w[k]*f(x[k])
    Gauss_values.append(s)
    
    x2,w2 = gaussxwab(2*N,a,b)
    for l in range(2*N):
        s2 += w2[l]*f(x2[l])
    Gauss_values_2N.append(s2)
    
list=[]

for j in range(3,12,1):
    list.append(j)
    print('The value for N=',2**j,'is:',Gauss_values[list.index(j)],'for Gaussian integration')
    
    print('                       ',simp_values[list.index(j)],'for Simpson integration')
    print('                       ',trap_values[list.index(j)],'for trapezoidal integration')
    print('')


simp_err = []
trap_err = []
Gauss_err_raw = []




for k in range(0,9,1):
    
    simp_err.append(abs((pi-simp_values[k]))/pi)
    trap_err.append(abs((pi-trap_values[k])/pi))
    Gauss_err_raw.append(abs((pi-Gauss_values[k])/pi))
    
   
plt.errorbar(N_values,simp_err,marker='o',label='relative error using Simpson')
plt.errorbar(N_values,trap_err,marker='o',label='relative error using trapezoidal')
plt.errorbar(N_values,Gauss_err_raw,marker='o',label='relative error using Gaussian')

plt.xscale('log')
plt.yscale('log')
plt.xlabel('Values of N (log scale)')
plt.ylabel('relative error (log scale)')
plt.legend()
plt.title('Relative errors using Simpson, trapezoidal and Gaussian integration')
plt.show()

Gauss_err_formula = []

for k in range(len(N_values)):
    Gauss_err_formula.append(abs((Gauss_values_2N[k] - Gauss_values[k])/pi))


plt.errorbar(N_values,Gauss_err_formula,marker='o',label='relative error using Gaussian')
plt.errorbar(N_values,trap_err,marker='o',label='relative error using trapezoidal')
plt.errorbar(N_values,simp_err,marker='o',label='relative error using Simpson')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Values of N (log scale)')
plt.ylabel('relative error (log scale)')
plt.legend()
plt.title('Relative errors using Simpson, trapezoidal and Gaussian integration')
plt.show()