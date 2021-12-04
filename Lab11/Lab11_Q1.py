# -*- coding: utf-8 -*-
"""
Created on Mon Nov 29 03:37:28 2021

@author: farzh
"""

from random import random,randrange
from math import exp,pi
from numpy import ones
from pylab import plot,ylabel,show
import numpy as np
import matplotlib.pyplot as plt

def Energy(T,steps):

 N = 1000
 

 # Create a 2D array to store the quantum numbers
 n = ones([N,3],int)

 # Main loop
 eplot = []
 E = 3*N*pi*pi/2
 for k in range(steps):

    # Choose the particle and the move
    i = randrange(N)
    j = randrange(3)
    if random()<0.5:
        dn = 1
        dE = (2*n[i,j]+1)*pi*pi/2
    else:
        dn = -1
        dE = (-2*n[i,j]+1)*pi*pi/2

    # Decide whether to accept the move
    if n[i,j]>1 or dn==1:
        if random()<exp(-dE/T):
            n[i,j] += dn
            E += dE

    eplot.append(E)

# # Make the graph
#  plot(eplot)
#  ylabel("Energy")
#  show()

# # This calculates the energy of each particle, neglecting constant factors
 energy_n = n[:, 0]**2 + n[:, 1]**2 + n[:, 2]**2
#  # This calculates the frequency distribution and creates a plot
#  plt.figure(2)
#  plt.clf()
 hist_output = plt.hist(energy_n, 50)
#  # This is the frequency distribution
#  energy_frequency = hist_output[0]
#  # This is what the x-axis of the plot should look like
#  # if we plot the energy distribution as a function of n
#  # the 2nd axis of hist_output contains the boundaries of the array.
#  # Instead, we want their central value, below.
 energy_vals = 0.5*(hist_output[1][:-1] + hist_output[1][1:])
 n_vals = energy_vals**0.5
#  # Create the desired plot
#  plt.figure(3)
#  plt.clf()
#  plt.bar(n_vals, energy_frequency, width=0.1)

 n_avg = np.sum(hist_output[0]*n_vals)/np.sum(hist_output[0])

 return n_avg, E

temps= (10,40,100,400,1200,1600)
step=(250000,250000,250000,650000,2000000,2500000)
Energies = []
n_avg =[]
for i in range(len(temps)):
    Energies.append(Energy(temps[i],step[i])[1])
    n_avg.append(Energy(temps[i],step[i])[0])
    
#print(Energies)  
#print(n_avg)  
plt.plot(temps,Energies,marker='o',color='blue')

plt.xlabel('Temperature(Kelvin)')
plt.ylabel('Energy')

plt.title('Energy vs temperature plot for non-interacting quantum ideal gas in a box')
plt.show()

plt.errorbar(temps,n_avg,marker='o')
plt.xlabel('Temperature(Kelvin)')
plt.ylabel('average "n" value')
plt.title('n vs temperature plot for non-interacting quantum ideal gas in a box')
plt.show()

Cv = (Energies[-1] - Energies[0])/(temps[-1] - temps[0])
print ("Over this range, heat capacity is %.2f" % Cv)
