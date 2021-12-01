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

from numpy.lib.function_base import append

N = 25
R = 0.02
Tmax = 10.0
Tmin = 1e-3
tau = 1e5
n = 5
compute_paths = True

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

# Set up the graphics
# display(center=[0.5,0.5])
# for i in range(N):
#     sphere(pos=r[i],radius=R)
# l = curve(pos=r,radius=R/2)


Ds = []
rs = []

if compute_paths:
    # Main loop
    for i in range (n):
        ns = (time.time() * 1000)*(i+1)
        seed(ns) # randomize the seed for the following paths
        #print(ns)
        t = 0
        T = Tmax
        while T>Tmin:

            # Cooling
            t += 1
            T = Tmax*exp(-t/tau)

            # Update the visualization every 100 moves
            # if t%100==0:
            #     l.pos = r

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
offset = 0.004
colorlist = ['red','orange','limegreen','blue','blueviolet']
for i in range(n):
    plt.plot(rs[i][:, 0]+offset*(i-n/2), rs[i][:, 1]+offset*(i-n/2), '-', c=colorlist[i%5], label='D = {0:.2f}'.format(Ds[i]))  # plot paths

plt.title('Travelling Salesman Paths,\n tau = {0:.1f}'.format(tau))
plt.legend()
# plot monomers

# plt.xlim([N/3.0, 5.0*N/3.0])
# plt.ylim([N/3.0, 5.0*N/3.0])
plt.axis('equal')
# plt.xticks([])  # we just want to see the shape
# plt.yticks([])
plt.tight_layout()

dpi = 150
plt.savefig('Q3a_final_path_Tau{0:.1f}_N{1:d}_n{2:}_Tmax{3:.1f}.png'.format(tau, N, n, Tmax),
            dpi=dpi)

# print('Energy averaged over last half of simulations is: {0:.2f}'
#       .format(np.mean(energy_array[n//2:])))
print('Total distance: {:.2f} +- {:.2f}'.format(np.average(Ds), np.std(Ds)))
plt.show()