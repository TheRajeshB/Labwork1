# -*- coding: utf-8 -*-
"""
Created on Fri Oct 29 14:57:59 2021

@author: farzh
"""
import numpy as np
import matplotlib.pyplot as plt
m = 9.1094e-31 # Mass of electron in kg
hbar = 1.0546e-34 # Planck's constant over 2*pi
e = 1.6022e-19 #charge of an electron in Coulomb

N = 1000 #no. of steps used in RK4 method 
a = 10**(-11)
 
x_start = -10*a #start of range over which wave-function is calculated 
x_end = 10*a #end of range over which wave-function is calculated 

V_0 = 50*e #constant to be used in calculation of potential energy (in Joules)

h = (x_end - x_start)/N #interval for RK4 method 

def V(x):
    ''' function that takes in position 'x' and outputs
        potential energy associated with that position'''
        
    return  V_0*(x**2)/(a**2)

def f(r,x,E):
    ''' function defining ODE's for wavefunction(psi) 
        and phi(derivative of psi wih respect to x)'''    
    psi = r[0] #wavefunction 
    phi = r[1] #rate of change of wave-function with respect to position
    f_psi = phi
    f_phi = (2*m/hbar**2)*(V(x)-E)*psi 
    
    return np.array([f_psi,f_phi],float)
wave_values = []
x_values = []
def wavefunction(E):
        ''' function calculates values of a wavefunction associated with some energy E'''
        psi = 0.0
        phi= 1.0
        r = np.array([psi,phi],float)
        for x in np.arange(x_start,x_end,h):
          x_values.append(x)         
          k1 = h*f(r,x,E) 
          k2 = h*f(r+0.5*k1,x+0.5*h,E) 
          k3 = h*f(r+0.5*k2,x+0.5*h,E) 
          k4 =h*f(r+k3,x+h,E) 
          r += (k1+2*k2+2*k3+k4)/6
          wave_values.append(r[0])
        return r[0]

def secant_root(E1,E2):
  ''' Using secant root method to find energy to find energy'''
   
  psi2 = wavefunction(E1)
  
  target_accuracy = e/1000
  while abs(E1 - E2)> target_accuracy:
    psi1,psi2 = psi2, wavefunction(E2)
    E1, E2 = E2, E2-psi2*(E2-E1)/(psi2-psi1)
        
  return E2/e 

print('Energy of ground state is:',secant_root(0,e),'eV')
print('Energy of first excited state is:',secant_root(200*e, 400*e),'eV')
print('Energy of second excited state is:',secant_root(500*e,700*e),'eV')


