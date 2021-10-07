# -*- coding: utf-8 -*-
"""
Created on Wed Oct  6 15:46:49 2021

@author: farzh
"""

from math import exp 
import numpy as np
import matplotlib.pyplot as plt
def f(c):
    steps = 1
    def g(x):
        return 1 - exp(-c*x)
    
     

    def err(x_current,x_next):

        derivative  = c*exp(-c*x_current)
        return (x_current - x_next)/(1 - (1/derivative))
    
    x_current = 1.0
    x_next = g(x_current)
    while abs(err(x_current,x_next)) > 10**-6:
        x_current =x_next
        x_next = g(x_next)
        steps += 1
    
   
    
    return x_next
    
    
        
        
    
f_values=[]
        
c = np.arange(0.01,3.01,0.01)
for i in c:
    f_values.append(f(i))

plt.plot(c,f_values,label='result using relaxation method')
plt.xlabel('c values')
plt.ylabel('x values')
plt.title('Plot showing results of relaxation method for various c-values')
plt.legend()
plt.show()