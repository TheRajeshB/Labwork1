# -*- coding: utf-8 -*-
'''Created on 2021-12-01

@author: Gabriel Bailey
Adapted from salesman.py

This code will test the sensitivity to tau of the simulated annealing optimization for the travelling salesman problem. It will output a plot of some of the paths taken.
'''
import matplotlib.pyplot as plt
import numpy as np

from math import sqrt,exp
from numpy import empty
from random import random,randrange,seed
import time

N = 25 # Number of cities
n = 5 # Number of different routes to find
Tmax = 10.0
Tmin = 1e-3
tau = 1e4

compute_paths = True # Compute values rather than use saved ones

ns = 14
seed(ns)

# Function to calculate the magnitude of a vector
def mag(x):
    return sqrt(x[0]**2+x[1]**2)

# Function to calculate the total length of the tour
def distance():
    s = 0.0
    for i in range(N):
        s += mag(r[i+1]-r[i])
    return s

# Choose N city locations and calculate the initial distance
r = empty([N+1,2],float)
for i in range(N):
    r[i,0] = random()
    r[i,1] = random()
r[N] = r[0]
D = distance()

# Lists for storing each annealing processes data (will have length n)
Ds = []
rs = []

if compute_paths:
    # Main loop
    for i in range (n):
        # Change the seed for each path
        #ns = (time.time() * 1000)*(i+1) # Really randomize it
        ns = ns*i
        seed(ns) 
        #print(ns)
        t = 0
        T = Tmax
        while T>Tmin:

            # Cooling
            t += 1
            T = Tmax*exp(-t/tau)

            # Choose two cities to swap and make sure they are distinct
            i,j = randrange(1,N),randrange(1,N)
            while i==j:
                i,j = randrange(1,N),randrange(1,N)

            # Swap them and calculate the change in distance
            oldD = D
            r[i,0],r[j,0] = r[j,0],r[i,0]
            r[i,1],r[j,1] = r[j,1],r[i,1]
            D = distance()
            deltaD = D - oldD

            # If the move is rejected, swap them back again
            if random()>exp(-deltaD/T):
                r[i,0],r[j,0] = r[j,0],r[i,0]
                r[i,1],r[j,1] = r[j,1],r[i,1]
                D = oldD
        Ds.append(D)
        rs.append(r.copy())

    np.savez('Q3a_paths_{:.2f}_{}'.format(tau,n), rs=rs, Ds=Ds,)
else:
    npzfile = np.load('Q3a_paths_{:.2f}_{}.npz'.format(tau,n))
    rs = npzfile['rs']
    Ds = npzfile['Ds']


plt.figure()

plt.scatter(r[:, 0], r[:, 1], marker='o', c='r')  # plot Cities
plt.plot(r[0, 0], r[0, 1], marker='o', c='k', label='Start/End')  # plot Start

# Try to plot each of the different paths, with an offset in x,y that will allow you to see each one
# (Doesn't work for lines parallel to x=y)
offset = 0.006
colorlist = ['red','orange','limegreen','blue','blueviolet']
colorlist.reverse()
for i in range(n):
    plt.plot(rs[i][:, 0], rs[i][:, 1]+offset*(i-n/2), '-', c=colorlist[i%5], label='D = {0:.2f}'.format(Ds[i]))  # plot paths

plt.title('Travelling Salesman Paths,\n tau = {0:.1f}'.format(tau))
plt.legend()
plt.axis('equal')
plt.tight_layout()

dpi = 150
plt.savefig('Q3a_final_path_Tau{0:.1f}_N{1:d}_n{2:}_Tmax{3:.1f}.png'.format(tau, N, n, Tmax),
            dpi=dpi)

print('Total distance: {:.2f} +- {:.2f}'.format(np.average(Ds), np.std(Ds)))
plt.show()