# -*- coding: utf-8 -*-
"""
Created on Fri Dec  3 19:01:20 2021

@author: farzh
"""

# import modules
import numpy as np
from random import random, randrange
import matplotlib.pyplot as plt



value = random()
# define constants
kB = 1.0
T = 1.0 #changed to 2 or 3 as needed
J = 1.0
Beta = 1/kB*T

N = 100
steps = 10**5

def energyfunction(J_, dipoles):
    """ function to calculate energy """
    energy = -J_*np.sum(dipoles[0:-1]*dipoles[1:])
    return energy
def acceptance(delta_E):
    """ Function that takes energy change between final and initial state\
         and calculates acceptance probability; completed """
    #calculating probability of acceptance     
    if delta_E <=0:
        P_a = 1
    else:
        P_a = np.exp(-delta_E * Beta)
    
    if value < P_a:
        result = True
    else: 
        result = False
    
    
    return result  # result is True of False

# generate array of dipoles and initialize diagnostic quantities
dipoles = np.empty([N, N], dtype = int)  # empty matrix to store dipole values


#randoming original diploe oritentation
for i in range(0,N):
    for j in range(0,N):
       flip = randrange(0,2,1) #generates either 0 or 1, randomly 
       if flip==0:
           dipoles[i,j] = 1
       if flip==1:
           dipoles[i,j] = -1
       
    
energy = []  # empty list; to add to it, use energy.append(value)
magnet = []  # empty list; to add to it, use magnet.append(value)
time = [0] #list to store number of time step ( starting at 0 for intial state)
E = energyfunction(J, dipoles)
energy.append(E)
magnet.append(np.sum(dipoles))
dipole_state=np.array([])
np.append(dipole_state,dipoles)

print('Initial Energy is:',E)
print('Initial magentization is:', magnet)

for i in range(steps):
    a = randrange(0,N,1) #picks random integer value between 0 and N-1 for row value
    b = randrange(0,N,1)#picks random integer value between 0 and N-1 for column value
    picked = randrange(N)  # choose a victim
    dipoles[picked] *= -1  # propose to flip the victim
    Enew = energyfunction(J, dipoles)  # compute Energy of proposed new state

    delta_e = Enew-E #calculating change in energy
   
    accepted = acceptance(delta_e)
    #if the flip is accepted, enacting said flip
    if accepted == True:
        if dipoles[a,b]==1: #if orinial orientation was 1 \
            dipoles[a,b] =-1 #flip it to -1
        if dipoles[a,b]==-1: #if original orientation was -1 \
            dipoles[a,b]=1#flip it to 1
        
    # store energy and magnetization
    energy.append(Enew)
    magnet.append(np.sum(dipoles))
    time.append(i)
    #generating animation at 100 step intervals
    if i%100==0:
        plt.clf()
        plt.imshow(dipoles)
        plt.title('Plot of spins of dipoles for T={} for {} steps'.format(T,i))
        plt.colorbar()
        plt.pause(0.01)
        
    


# plot energy, magnetization

# plt.plot(time,energy)
# plt.title('Energy vs number of steps for 10^5 steps, at T={}'.format(T))
# plt.xlabel('Number of steps')
# plt.ylabel('Energy')
# plt.show()

# plt.plot(time,magnet)
# plt.title('Energy vs number of steps for 10^5 steps, at T={}'.format(T))
# plt.xlabel('Number of steps')
# plt.ylabel('Magnetization value')
# plt.show()

plt.imshow(dipoles)
plt.title('Plot of spins of dipoles for T={} for 10^5 steps(final state)'.format(T))
plt.colorbar()
plt.show()


 