# -*- coding: utf-8 -*-
'''Created on 2021-11-23

@author: Gabriel Bailey

This code will simulate photons leaving the photosphere of a star to generate the angle at which they leave and the limb-darkening effect. It will produce 2 plot: One for the final scattering angle and one for the intensity. The variable 'part' determines whether the results for part b or part c are generated.
'''
from logging import error
from random import randrange, uniform
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

# The function we are integrating over
def func(x):
    return(x**(-1/2)/(1+np.exp(x)))

# The weight function for importance sampling
def w(x):
    return(x**(-1/2))

# f/w
def f_over_w(x):
    return(1/(1+np.exp(x)))

# Probability distrobution used
def p(x):
    return 1/2 * x**(-1/2)

# Calculate the integral of f from a to b using mean value monte carlo with N iterations
def integrate_mv(f,a,b,N):
    f_sum = 0
    for i in range(N):
        x = (b-a)*np.random.random()
        #print(x)
        f_sum += f(x)
    #print(f_sum,a,b,N)
    #print((b-a)/N * f_sum)
    return (b-a)/N * f_sum

# Calculate the integral of f from a to b using importance sampling monte carlo with N iterations
def integrate_is(f,a,b,N):
    f_sum = 0
    for i in range(N):
        x = (np.random.random())**2
        
        f_sum += f(x)
    # print(f_sum,N)
    # print(1/N * f_sum)
    return 2/N * f_sum

#Constants


# options
compute_mv = True
compute_is = True
part = 'a' # Choose between part b and c
N = 10000
repeats = 100

f_mv = func
f_is = f_over_w
a = 0
b = 1
if part == 'a':
    f_mv = func
    f_is = f_over_w
    a = 0
    b = 1    
elif part == 'b':
    f_mv = func
    f_is = f_over_w
else:
    error('Invalid part of question selected (must be \'a\' or \'b\').')

# Array to put data into
I_mv = np.empty(repeats)

I_is = np.empty(repeats)

# Save data rather than recalculate it if desired:
if compute_mv:
    for i in range(repeats):
        I_mv[i] = integrate_mv(f_mv,a,b,N)
        #print(I_mv[i])

    np.savez('Q3_mv_'+part, I_mv=I_mv)
else:
    npzfile = np.load('Q3_mv_'+part+'.npz')
    I_mv = npzfile['I_mv']

if compute_is:
    for i in range(repeats):
        I_is[i] = integrate_is(f_is,a,b,N)

    np.savez('Q3_is_'+part, I_is=I_is)
else:
    npzfile = np.load('Q3_is_'+part+'.npz')
    I_is = npzfile['I_is']

# print('I_mv',I_mv)
# print('I_is',I_is)
# Plot the resulting mu histogram
fig, axmv = plt.subplots()
axmv.set_title('Mean Value Integral Spread') 
axmv.set_xlabel('Calculated Value')
axmv.set_ylabel('N') 
axmv.hist(I_mv, 10, range=[0.8, 0.88])

fig, axis = plt.subplots()
axis.set_title('Importance Sampling Integral Spread') 
axis.set_xlabel('Calculated Value')
axis.set_ylabel('N') 
axis.hist(I_is, 10, range=[0.8, 0.88])
# x = np.linspace(0,1,100)
# y = func(x)
# plt.plot(x,y)

# x = np.linspace(0,1,100)
# y = f_over_w(x)
# plt.plot(x,y)

plt.show()