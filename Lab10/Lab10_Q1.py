# -*- coding: utf-8 -*-
"""
Created on Thu Nov 25 22:20:57 2021

@author: farzh
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import RegularGridInterpolator
#loading data
loaded = np.load('Earth.npz')
data = loaded['data'] #array for land/ocean points
lon_array = loaded['lon'] #array for longitude values
lat_array = loaded['lat'] #array for latitude values


#The followng is the solution to Q1b

N= 5000 #number of points 

#generating random values for x and y 
a= np.random.random(N)
b = np.random.random(N)


#finding theta and phi from formulae

theta = np.arccos(1-2*a) 
phi = 2*np.pi*b


#defining Cartesian co-ordinates on globe 
x = np.sin(theta)*np.cos(phi)
y = np.sin(theta)*np.sin(phi)
z = np.cos(theta)

#generating plot 
fig = plt.figure()
ax = plt.axes(projection ='3d')
ax.scatter(x,y,z,s=1)
ax.set_xlabel('X-axis')
ax.set_ylabel('Y-axis')
ax.set_zlabel('Z-axis')
plt.title('Plot of randomly generated theta and phi points')
plt.show()

#The followng is the solution to Q1c
total_land = 0

for i in data: 
    'loop to find total number of land points'
    z= np.sum(i) #finding total number of land points in each row of matrix
    total_land += z #adding up total number of land points over all rows
    
#finding total number of points 
total = np.shape(data)[0]*np.shape(data)[1] #total number of points (land + sea)
land_frac = total_land/total #  land fraction
print('The land fraction for Q1c is :',land_frac)


#The followng is the solution to Q1d


lon_array = np.append(lon_array, 180) #adding value of 180 to longitude array to ensure all calculated(randomised) coordinates fall within bounds of this array 
c = np.copy(data[-1]) #copying last row of values from land/water array
d = np.array([c])
grid=np.vstack((data,d)) #setting the values for land/water along 180 longitude as the same as last rowof original 'data' array 

def Land_fraction(N):
    'function that takes in number of points(N),\
     generates N random points, finds the associated latitude/longtiude\
     determines if the point is on land or water and returns the fraction of land points'
     
    
    a= np.random.random(N) #generating random points
    b = np.random.random(N) #generating random points
     
    theta1 = np.arccos(1-2*a) - (np.pi/2) #generating theta with range same as that of latitude
    phi1 = 2*np.pi*b - np.pi #generating phi with range same as that of longitude 

    lat_random = theta1*(180/np.pi)  #converting to degrees to give randomly generated latitude 
    long_random = phi1*(180/np.pi)  #converting to degrees to give randomly generated longitude
    
    #settin up conditions for nearest interpolation
    
    interp = RegularGridInterpolator((lon_array, lat_array), grid, method='nearest')
    
    #lists to store values
    array=[]

    # land_theta = []
    land_phi=[]
    # water_theta = []
    # water_phi = []
    land_i = []
    water_i = []
    for i in range(N):
     'this loop carries out nearest interpolation for all randomly generated points'
     x1,y1 = long_random[i],lat_random[i]
     nearest_value = interp([x1, y1])
     array.append(nearest_value)
     
     
     if 1-10**-5<=nearest_value <= 1+10**-5:
         land_i.append(i)
         land_phi.append(y)
     elif -10**-5<=nearest_value <=10**-5: 
          water_i.append(i)
          

     
    
    z= 0 
    for i in array:
     z+= i #Finding total number of random points that are land points
    
    return z/N ,land_i,water_i,lat_random,long_random

print('The calculated land fraction for N=50 is:', Land_fraction(50)[0])
print('The calculated land fraction for N=500 is:', Land_fraction(500)[0])
print('The calculated land fraction for N=5000 is:', Land_fraction(5000)[0])
print('The calculated land fraction for N=50000 is:', Land_fraction(50000)[0])




land_i = Land_fraction(50000)[1]
water_i = Land_fraction(50000)[2]
lat_random = Land_fraction(50000)[3]
long_random = Land_fraction(50000)[4]
lat_land= []
long_land = []
lat_water=[]
long_water=[]
for i in land_i:
    lat_land.append(lat_random[i])
    long_land.append(long_random[i])
for i in water_i :
    lat_water.append(lat_random[i])
    long_water.append(long_random[i])

for i in range(len(long_land)):
  plt.plot(long_land[i],lat_land[i],marker='o',color = 'brown',markersize=0.5)
  plt.plot(long_water[i],lat_water[i],marker='o',color = 'blue',markersize=0.2)

plt.xlabel('longitude(degrees)')
plt.ylabel('latitude(degrees)')
plt.title('Land/water plot')

plt.show()
