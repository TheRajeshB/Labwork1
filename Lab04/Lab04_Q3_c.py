# -*- coding: utf-8 -*-
"""
Created on Wed Oct  6 20:34:55 2021

@author: farzh
"""
from math import exp
import numpy as np
from scipy.constants import h,k,c

def binary_method(x1,x2):
    """function that takes set of points and uses binary method to 
       calculate value of root for function
       x1 = left point of bracket around root
       x2 = right end of bracket around root """
    
    def f(x):
        """funcion for which root is to be calculated """
        
        return 5*(exp(-x))+x   - 5
    accuracy = 10**-6 #desired accuracy for value of root

    
    while abs(x1 -x2)> accuracy:
      if f(x1)==0: #checking if x1 is root 
          x2=x1
          x=x1 #re-setting values so as to end calculation
      if f(x2)==0:#checking if x2 is root
          x1 =x2
          x=x2 #re-setting values so as to end calculation
      
       
      if f(x1)>0 and f(x2)<0 or f(x1)<0 and f(x2)>0 :#check if sign of value of function at ends of bracket are different
        x = (x1+x2)/2 #calculating mid-point
      
        if np.sign(f(x)) == np.sign(f(x1)): #checking if value of function at mid-point is same as at x1
          x1 = x #setting new value for x1
      
        elif np.sign(f(x)) == np.sign(f(x2)): #checking if value of function at mid-point is same as at x2
          x2 = x #setting new value for x2
      

      elif np.sign(f(x1))==np.sign(f(x2)):#printing warning to adjust imput points if function has same sign at both points
          print('function output has same sign for both points')
          
      
          
      
    return (x1+x2)/2
          
    
x = binary_method(1,10)
b = h*c/(k*x) #calculating value of displacement constant
print('displacement constant is:',b,'mK')
lamb = 502*(10**-9) #wavelength is meter
print('surface temperature of the sun is:',b/lamb,'Kelvin')
