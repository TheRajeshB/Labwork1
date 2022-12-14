# -*- coding: utf-8 -*-
'''Created on 2021-11-23

@author: Gabriel Bailey

This code will integrate x^(-1/2)/(1+e^x) from 0 to 1 using monte carlo methods, specifically mean-value and importance sampling. It will produce 2 histograms of the result, one for each method.
'''
import matplotlib.pyplot as plt
import numpy as np

# The function we are integrating over
def f(x):
    return(x**(-1/2)/(1+np.exp(x)))

# f/w
def f_over_w(x):
    return(1/(1+np.exp(x)))

# Calculate the integral of f from a to b using mean value monte carlo with N iterations
def integrate_mv(f,a,b,N):
    f_sum = 0
    for i in range(N):
        x = (b-a)*np.random.random()
        f_sum += f(x)
    return (b-a)/N * f_sum

# Calculate the integral of f from a to b using importance sampling monte carlo with N iterations
def integrate_is(f,a,b,N):
    f_sum = 0
    for i in range(N):
        # Generate x based on the probability distrobution p
        x = (np.random.random())**2
        
        f_sum += f(x)
    # The 2 is the integral of w from a to b
    return 2 * 1/N * f_sum

# Options
compute_mv = True
compute_is = True

# Constants
N = 10000
repeats = 100

f_mv = f
f_is = f_over_w
a = 0
b = 1

# Arrays to put data into
I_mv = np.empty(repeats)
I_is = np.empty(repeats)

# Save data rather than recalculate it if desired:
if compute_mv:
    for i in range(repeats):
        I_mv[i] = integrate_mv(f_mv,a,b,N)

    np.savez('Q3_mv_a', I_mv=I_mv)
else:
    npzfile = np.load('Q3_mv_a.npz')
    I_mv = npzfile['I_mv']

if compute_is:
    for i in range(repeats):
        I_is[i] = integrate_is(f_is,a,b,N)

    np.savez('Q3_is_a', I_is=I_is)
else:
    npzfile = np.load('Q3_is_a.npz')
    I_is = npzfile['I_is']

# Plot the resulting mu histogram
fig, axmv = plt.subplots()
axmv.set_title('Mean Value Integral Histogram') 
axmv.set_xlabel('Calculated Value')
axmv.set_ylabel('Count') 
axmv.hist(I_mv, 10, range=[0.8, 0.88])

fig, axis = plt.subplots()
axis.set_title('Importance Sampling Integral Histogram') 
axis.set_xlabel('Calculated Value')
axis.set_ylabel('Count') 
axis.hist(I_is, 10, range=[0.8, 0.88])

plt.show()