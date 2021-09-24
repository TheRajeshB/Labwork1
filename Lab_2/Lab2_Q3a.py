# -*- coding: utf-8 -*-
"""
Created on Tue Sep 21 13:21:29 2021

@author: farzh
"""
from math import *
import numpy as np
from scipy.integrate import quad
from scipy.special import kn
from IntegrationFunctions import symb_integrate, trap_integrate, simp_integrate
import matplotlib.pyplot as plt
#defining various terms and functions used in eqn 8 



def V(u):
    Q =10**-13 #charge in Couloumb
    m = tan(u) #defined to simplify formula
    n = cos(u) #defined to simplify formula
    z = 0 #dist along z-axis (z) in m
    epsilon = 8.854e-12 #epsilon-nought in calculation of electrostatic constant
    l = 0.001 # in mm
    k = Q/(4*pi*epsilon)
   
    
    r = np.linspace(0.00025,0.005,100)  #range of values for 'r' in mm
    c = ((z-l*m)**2 + r**2)**0.5 #terms under root in denominator
    
    
    return (k*e**-(m**2))/((n**2)*c)
    

#fix to simpsons
def simp(f,a,b,N):
    '''function to perform Simspon integration
    a = start point of integration
    b = end point of integration
    N = no.of slices'''
    
  
    h= (b-a)/N #length of each interval
    
    sum_odd = 0.0
    sum_even = 0.0
    initial_sum = f(a) + f(b) 
    
    trap_area=0.0
    
    #calculating last slice in trapezoidal if N is odd 
    if N % 2 == 1: #if N is odd calc the last slice with trapezoidal rule
        trap_area += h*f(a+(N-1)*h)
        N -= 1
        b=b-h #setting new end point of Simpson integration to exclude area covered by trapezoidal integration
    
    
    
    #loop to find sum of odd terms for Simpson integration
    for k in range(1,int(N),2):
        sum_odd  += f((a+k*h))
     
    #loop to find sum of even terms for Simpson integration
    for k in range(2,int(N),2):
        sum_even += f((a+k*h))
   
    return(h*( initial_sum + 4*sum_odd + 2*sum_even)/3) + (trap_area)

#assigning values to terms to be used in eqn 9
Q =10**-13 #charge in Couloumb
epsilon = 8.854e-12 #epsilon-nought in calculation of electrostatic constant
l=0.001 #in m
r = np.linspace(0.00025,0.005,100)  #range of values for 'r' in m  
k = Q/(4*pi*epsilon) 

d= (r**2)/(2*(l**2)) #input for Bessel function
#calculating eqn 9 for range of 'r' values
eqn_9 = k*(e**d)*kn(0,d)/l



plt.errorbar(r,simp(V,-pi/2,pi/2,10),label='Simpson integration result for N=10')
plt.errorbar(r,eqn_9,label='Result from equation 9(overlaps for N>8)')
plt.xlabel('radial distance (=r in mm)')
plt.ylabel('potential at r (Volts)')
plt.title('graph showing comparison of results between Simpson integration and eqn 9')
plt.legend()
plt.show()

#first value for fraction error less than 1 in a million is 53
diff = eqn_9 -simp(V,-pi/2,pi/2,10) 
plt.plot(r,diff,label='Difference for N =10')
plt.axhline(y=0, xmin=0.0, xmax=5.0, color='r', linestyle='-.', linewidth=3,label='0-error line')
plt.xlabel('radial distance(=r in mm)')
plt.ylabel('Difference ')
plt.legend()
plt.title('graph showing difference between results from eqn 9 and Simpson integration')
plt.show()


#checking fractional error:
for i in (diff/eqn_9) :
    if i >10**-6:
     print('fractional error is over 1 in a million')
    else:
     None
    