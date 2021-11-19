# -*- coding: utf-8 -*-
"""
Created on Fri Nov 19 11:24:04 2021

@author: farzh
"""

import numpy as np 
import matplotlib.pyplot as plt

#setting constants
epsilon = 1
delta_t = 0.005
delta_x = 0.02
Lx = 2*np.pi
Tf = 2

Nx = int(Lx//delta_x)
Nt = int(Tf//delta_t)

#matrix to hold values 
u = np.zeros([Nx, Nt],float)

# Set boundary conditions
u[0,:] = 0
u[Nx-1,:] = 0
#set initial conditions
for x in range(1,Nx-1):
    u[x,0] = np.sin(x*delta_x)
    # Also forward Euler step:
    u[x,1] = u[x,0] - epsilon*delta_t/(4*delta_x)*(u[x+1,0]**2-u[x-1,0]**2)    
beta = epsilon*delta_t/delta_x

for t in range(1,Nt-1):
    for x in range(1,Nx-1):
        a = ((u[x+1,t]**2) - (u[x,t]**2))
        b = ((u[x-1,t]**2) - (u[x,t]**2))
        u[x,t+1] = u[x,t] - beta*((u[x+1,t]**2) - (u[x-1,t]**2))/4 + \
            (beta**2)*((u[x,t] + u[x+1,t])*a + (u[x,t] + u[x-1,t])*b)/8 #implenting Lax-Wendroff method

def plot_result(u,time):
    t = int(time/Tf * Nt)
    #print(t)
    fig = plt.figure(figsize=[10,5])
    ax = fig.add_subplot(1,1,1)
    ax.set_xlim(0, Lx)
    ax.set_ylim(-1.6, 1.6)
    x = np.linspace(0, Lx, Nx)
    ax.plot(x,u[:,t])
        
    ax.grid(which='both', axis='y')
    ax.set_title('Burger\'s Equation at Time {}s using Lax-Wendroff'.format(time)) 
    ax.set_xlabel('x')
    ax.set_ylabel('u')
    return ax
plot_result(u,0.0)
plot_result(u,0.5)
plot_result(u,1.0)
plot_result(u,1.5)
