# -*- coding: utf-8 -*-
'''Created on 2021-12-02

@author: Gabriel Bailey
Adapted from Q3a

This code will find the minium value of a function via annealing. It will print data about the simulation including it's final values, and then output plots of the x, y, and f(x,y) values over time.
'''
import matplotlib.pyplot as plt
import numpy as np

from math import sqrt,exp,cos,pi
from random import random,seed

# For i : 10.0, 1e-4, 1e5
# For ii: 10.0, 1e-4, 1e5 (same)
Tmax = 10.0
Tmin = 1e-4
tau = 1e5
compute_data = True
part = 'ii'

ns = 20
seed(ns)

# Function to generate random gaussian number in a repeatable way
def gauss_dist(sigma):
    z = random()
    theta = 2*pi*random()
    r = (-2*sigma*np.log(1-z))**(1/2)
    x = r*np.cos(theta)
    y = r*np.sin(theta)
    return x,y

# The function we are finding the min of
def func(x,y):
    if part == 'i':
        return x**2 - cos(4*pi*x) + (y-1)**2
    elif part == 'ii':
        return cos(x) + cos(sqrt(2)*x) + cos(sqrt(3)*x) + (y-1)**2
    else:
        raise ValueError('Wrong part: only parts i and ii supported')

t = 0
xs = [2]
ys = [2]
f_vals = [func(xs[0], ys[0])]

if compute_data:
    # Main loop
    
    T = Tmax
    while T>Tmin:

        # Cooling
        t += 1
        T = Tmax*exp(-t/tau)

        # Monte Carlo move
        dx,dy = gauss_dist(1)
        x = xs[t-1]+dx
        y = ys[t-1]+dy
        f_val = func(x,y)
        
        delta_f = f_val - f_vals[t-1]
        
        # Reject values outside the given range for part ii
        extra_cond = (part == 'ii') and not ((0<x<50) and (-20<y<20)) 
        # If the move is rejected, reuse the old values, otherwise add the new ones
        if random()>exp(-delta_f/T) or extra_cond:
            xs.append(xs[t-1])
            ys.append(ys[t-1])
            f_vals.append(f_vals[t-1])
        else:
            xs.append(x)
            ys.append(y)
            f_vals.append(f_val)

    # Save or load depending on whether compute_data is set to true or not
    np.savez('Q3b{}_data_{:.2f}_{:.2f}_{:.2f}'.format(part,tau,Tmax,Tmin), xs=xs, ys=ys, f_vals=f_vals, t=t)
else:
    npzfile = np.load('Q3b{}_data_{:.2f}_{:.2f}_{:.2f}.npz'.format(part,tau,Tmax,Tmin))
    xs = npzfile['xs']
    ys = npzfile['ys']
    f_vals = npzfile['f_vals']
    t = npzfile['t']

# Print the simulation info
print('Paramaters: Tau={}, Tmax={}, Tmin={}'.format(tau,Tmax,Tmin))

print('Final f value: {:}'.format(f_vals[-1]))
print('max f value: {:}'.format(np.max(f_vals)))
print('min f value: {:}'.format(np.min(f_vals)))

print('Final x y: {:f}, {:f}'.format(xs[-1],ys[-1]))

# Plot the results
dpi = 150
# Plot the f value over time (for testing)
""" plt.figure()
plt.title('Annealing Function Values,\n tau = {0:.1f}'.format(tau))
plt.xlabel('Time (steps)')
plt.ylabel('Function Value')
plt.plot(np.arange(len(f_vals)),f_vals, '.', )#label='D = {0:.2f}'.format(Ds[i]))  # plot paths
plt.tight_layout() """

# Plot the x value over time
plt.figure()
plt.title('Annealing x Values')
plt.xlabel('Time (steps)')
plt.ylabel('x Value')
plt.plot(np.arange(len(xs)),xs, ',', )
plt.savefig('Q3b{}_x.png'.format(part), dpi=dpi)

# Plot the y value over time
plt.figure()
plt.title('Annealing y Values')
plt.xlabel('Time (steps)')
plt.ylabel('y Value')
plt.plot(np.arange(len(ys)),ys, ',', )
plt.savefig('Q3b{}_y.png'.format(part), dpi=dpi)

plt.show()