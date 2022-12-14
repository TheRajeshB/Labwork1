# -*- coding: utf-8 -*-
"""
Created on Tue Nov  2 11:25:35 2021

@author: farzh
"""

import numpy as np
import matplotlib.pyplot as plt

#defining boundary 
a = 0.1 #grid spacing in cm 
L = int(20/a) #number of grid-points for length of domain 
W = int(8/a)  #number of grid points for width of domain  

#length of sides of domain 
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

phi = np.zeros([W+1,L+1],float)

delta = 1.0 
target = 10**-6 #target accuracy 




# Initialize boundary conditions
for i in range(W+1):
    for j in range(L+1):
        if i==0 and 0<=j<=L:
            phi[i,j]= 10 #Boundary condition for HG
        
        elif j==0 and 0<i<W:
            phi[i,j] = BC_HA[i] #Boundary condition for HA
            
        elif j==200 and 0<i<W:
            phi[i,j] = BC_GF[i] #Boundary condition for GF
            
        elif i == 80 and 0<j<50:
            phi[i,j] = BC_AB[j] #Boundary condition for AB
            
        elif j == 50 and 50<i<80:
            phi[i,j] = BC_CB[i-50] #Boundary condition for CB
            
        elif i == 50 and 50<j<150:
            phi[i,j] = 7 #Boundary condition for CD
        
        elif j==150 and 50<i<80:
            phi[i,j] = BC_DE[i-50] #Boundary condition for DE
            
        elif i==80 and 150<j<200:
            phi[i,j] = BC_EF[j-150] #Boundary condition for EF


def results(iterations,omega):
 """functions that takes in number of iterations and omega, performs Gauss-Siedel relaxation with
      replacements and over-realaxation for 
     that number of iterations, and outputs result (in matrix form)  """  
 for n in range(iterations):
    
    delta = 0.0
    for i in range(W):
        
        for j in range(L):
            if i==0 and 0<=j<=L:
                phi[i,j]= 10 #Boundary condition for HG
            
            elif j==0 and 0<i<W:
                phi[i,j] = BC_HA[i] #Boundary condition for HA
                
            elif j==200 and 0<i<W: 
                phi[i,j] = BC_GF[i] #Boundary condition for GF
                
            elif i == 80 and 0<j<50:
                phi[i,j] = BC_AB[j] #Boundary condition for AB
                
            elif j == 50 and 50<i<80:
                phi[i,j] = BC_CB[i-50] #Boundary condition for CB
                
            elif i == 50 and 50<j<150:
                phi[i,j] = 7 #Boundary condition for CD
            
            elif j==150 and 50<i<80:
                phi[i,j] = BC_DE[i-50] #Boundary condition for DE
                
            elif i==80 and 150<j<200:
                phi[i,j] = BC_EF[j-150] #Boundary condition for EF
                
            else:
                phi[i, j] =  (1 + omega) * (phi[i+1, j] + phi[i-1, j] + phi[i, j+1] + phi[i, j-1])/4 - omega * phi[i, j]
    
 
 return np.flip(phi,0) #adjusting matrix to make origin at bottom left of domain

plt.pcolormesh(results(100,0.9),cmap='Greens')
plt.title('Temperature distribution for omega =0.9')
plt.xlabel('x(mm)')
plt.ylabel('y(mm)')
plt.colorbar()
plt.show( )

for x in np.arange(0,100,9):
    plt.clf()
    plt.pcolormesh(results(x,0.9),cmap='Greens')
    plt.title('Tempreature distribution')
    plt.xlabel('x(mm)')
    plt.ylabel('y(mm)')
    plt.colorbar()
    plt.draw()
    plt.pause(0.1)
    
