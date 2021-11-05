# -*- coding: utf-8 -*-
"""
Created on Tue Nov  2 11:25:35 2021

@author: farzh
"""

import numpy as np
import matplotlib.pyplot as plt

#defining boundary 
a = 0.1 #grid spacing in cm 
l = 20
w = 8
L = int(l/a) #number of grid-points for length of domain 
W = int(w/a)  #number of grid points for width of domain  

#length of sides of shape 
AB = 5 #in cm
BC = 3 #in cm 
CD = 10 #in cm 
DE = 3 # in cm 
EF = 5 #in cm 
FG = 8 #in cm 
GH = 20 # in cm 
HA = 8 #in cm 

#defining boundary conditions 
BC_AB = np.linspace(0,5,int(AB/a)+1)
BC_CB = np.linspace(7,5,int(BC/a)+1)
BC_DE = np.linspace(7,5,int(DE/a)+1)
BC_EF = np.linspace(5,0,int(EF/a)+1)
BC_GF = np.linspace(10,0,int(FG/a)+1)
BC_HA = np.linspace(10,0,int(HA/a)+1)

N = 100 # number of iterations

phi = np.zeros([N,W+1,L+1],float)

delta = 1.0 
target = 10**-6

omega = 0.9
print(phi.shape)
phi[0,0,100]
# Initialize boundary conditions
for i in range(W+1):
    for j in range(L+1):
        if i==0 and 0<=j<=L:
            phi[0,i,j]= 10 #Boundary condition for HG
        
        elif j==0 and 0<=i<=W:
            phi[0,i,j] = BC_HA[i] #Boundary condition for HA
            
        elif j==L and 0<=i<=W:
            phi[0,i,j] = BC_GF[i] #Boundary condition for GF
            
        elif i == 80 and 0<=j<=50:
            phi[0,i,j] = BC_AB[j] #Boundary condition for AB
            
        elif j == 50 and 50<=i<=W:
            phi[0,i,j] = BC_CB[i-50] #Boundary condition for CB
            
        elif i == 50 and 50<=j<=150:
            phi[0,i,j] = 7 #Boundary condition for CD
        
        elif 50<=j<=150 and 50<=i<=W:
            phi[0,i,j] = BC_DE[i-50] #Boundary condition for DE
            
        elif i==80 and 150<=j<=L:
            phi[0,i,j] = BC_EF[j-150] #Boundary condition for EF
        else:
            #print(0,i,j)
            phi[0,i, j] =  (1 + omega) * (phi[0,i+1, j] + phi[0,i-1, j] + phi[0,i, j+1] + phi[0,i, j-1])/4 - omega * phi[0,i, j]
print(np.max(phi[0]),np.min(phi[0]))
# Run simulation

for n in range(1,N):
    delta = 0.0
    for i in range(W):
        for j in range(L):
            phi[n] = phi[n-1]
            if i==0 and 0<j<L:
                phi[n,i,j]= 10 #Boundary condition for HG
            
            elif j==0 and 0<i<W:
                phi[n,i,j] = BC_HA[i] #Boundary condition for HA
                
            elif j==200 and 0<i<W: 
                phi[n,i,j] = BC_GF[i] #Boundary condition for GF
                
            elif i == 80 and 0<j<50:
                phi[n,i,j] = BC_AB[j] #Boundary condition for AB
                
            elif j == 50 and 50<i<80:
                phi[n,i,j] = BC_CB[i-50] #Boundary condition for CB
                
            elif i == 50 and 50<j<150:
                phi[n,i,j] = 7 #Boundary condition for CD
            
            elif j==150 and 50<i<80:
                phi[n,i,j] = BC_DE[i-50] #Boundary condition for DE
                
            elif i==80 and 150<j<200:
                phi[n,i,j] = BC_EF[j-150] #Boundary condition for EF
                
            else:
                phi[n,i, j] =  (1 + omega) * (phi[n,i+1, j] + phi[n,i-1, j] + phi[n,i, j+1] + phi[n,i, j-1])/4 - omega * phi[n,i, j]

phi = np.flip(phi,0)

# fig = plt.figure(figsize=[10,5])
# ax = fig.add_subplot(1,1,1)
# #im = ax.imshow(m1, interpolation='None')
# x_coords = np.linspace(0,20,L+1)
# y_coords = np.linspace(0,8,W+1)
# ax.contourf(x_coords,y_coords,phi[0])
# ax.set_title('Temperature Distrobution, omega = {}'.format(omega)) 
# ax.set_xlabel('x (cm)')
# ax.set_ylabel('y (cm)')
# #plt.colorbar(ax)
# plt.show()

from matplotlib.animation import FuncAnimation, FFMpegWriter

fig = plt.figure()
ax = plt.axes(xlim=(0, l), ylim=(0, w), xlabel='x', ylabel='y')
#im = ax.imshow(m1, interpolation='None')
x = np.linspace(0,20,L+1)
y = np.linspace(0,8,W+1)
x,y = np.meshgrid(x,y)
cvals = np.linspace(0,10,N+1) 
cont = plt.contourf(x,y,phi[0],cvals)
plt.colorbar()

plt.title('Temperature Distrobution, omega = {}'.format(omega)) 
plt.xlabel('x (cm)')
plt.ylabel('y (cm)')
time_text = plt.text(0.1, 0.1, 'time = ', transform=ax.transAxes)

def update(i):
    global cont,x,y
    #print(x,y)
    for c in cont.collections:
        c.remove()  # removes only the contours, leaves the rest intact
    cont = plt.contourf(x, y, phi[i,:,:],cvals)
    #time_text.set_text('frame = %i' % i)
    return cont
    i = int(i)
    x = linspace(0, Lx, Nx)
    y = u[:,i]
    plot.set_data(x, y)
    time_text.set_text('time = {:1.2f}'.format(i*delta_t))
    return plot, time_text

ani = FuncAnimation(fig, update, frames=N, blit=True)

# Save the animation (need ffmpeg)
# f = 'Q3_vid.mp4' 
# writervideo = FFMpegWriter(fps=60) 
# ani.save(f, writer=writervideo)

plt.show()