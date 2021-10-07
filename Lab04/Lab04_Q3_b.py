# -*- coding: utf-8 -*-
"""
Created on Wed Oct  6 16:58:32 2021

@author: farzh
"""

from math import exp 
import numpy as np


def relaxation(c):
    """function takes in value of 'c'
    Uses relalaxation method for appropriate function
    Returns x-value that solves function and 
    number of steps that were needed to get answer """
    
    def g(x):
        """" function describing equation to be solved """
        return 1 - exp(-c*x)
    
     
    steps = 1 #starting number of iterations
    x_current = 1.0 #starting estimate
    x_next = g(x_current) #next estimate
    while abs(x_current-x_next) > 10**-6: #setting up threshold
        x_current =x_next #setting values for next iteration
        x_next = g(x_current)
        steps += 1 #increasing count of interations by 1
    
   
    
    return x_next, steps

print('the result and number of iterations using relaxation method is:',relaxation(2))

def over_relax(w,c):
    """function takes in value of 'c'
    Uses over-relalaxation method for appropriate function
    Returns x-value that solves function and 
    number of steps that were needed to get answer """
    def g(x):
        """" function describing equation to be solved """
        return 1- exp(-c*x)
    

    
    
    steps = 1#starting number of iterations
    x_current = 1.0
    x_next = x_current + (1+w)*(g(x_current)- x_current)
    
    while abs((x_current-x_next)) > 10**-6:
        x_current =x_next #setting values for next iteration
        x_next = x_current + (1+w)*(g(x_current)- x_current)
        steps += 1   #increasing count of interations by 1
    return x_next, steps

print('the result and number of iterations using relaxation method is:',over_relax(0.85,2),',for omega = 0.9')        
