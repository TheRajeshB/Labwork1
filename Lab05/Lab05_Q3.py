# -*- coding: utf-8 -*-
"""
Created on Tue Oct 12 21:04:57 2021

@author: farzh
"""

import numpy as np
import matplotlib.pyplot as plt

#loadng in data and assigning variables to store data
SLP = np.loadtxt('SLP.txt')
Time = np.loadtxt('times.txt')
long = np.loadtxt('lon.txt')

#creating zero-matrix with same dimensions as SLP
SLP_3 = np.zeros(SLP.shape)


for i in range(SLP.shape[0]): 
    
    SLP3_fourier = np.fft.rfft(SLP[i,:]) #finding Fourier transform for every row in SLP
    
    for j in range(SLP3_fourier.shape[0]):
        if j != 3: 
            SLP3_fourier[j] = 0#setting every row in Fourier transform of SLP not
                    #corresponding to Fourier wave number equal to 3 as zero
     
        
    SLP3_fourier_inv = np.fft.irfft(SLP3_fourier) #taking inverse Fourier transform of result of above computation
    SLP_3[i,:] = SLP3_fourier_inv #setting each output from above line as corresponding row 
                                  #in zero-matrix
    

SLP_5 = np.zeros(SLP.shape)


for i in range(SLP.shape[0]): 
    
    SLP5_fourier = np.fft.rfft(SLP[i,:]) #finding Fourier transform for every row in SLP
    
    for j in range(SLP5_fourier.shape[0]):
        if j != 5: 
            SLP5_fourier[j] = 0#setting every row in Fourier transform of SLP not
                    #corresponding to Fourier wave number equal to 5 as zero
     
        
    SLP5_fourier_inv = np.fft.irfft(SLP5_fourier) #taking inverse Fourier transform of result of above computation
    SLP_5[i,:] = SLP5_fourier_inv #setting each output from above line as corresponding row 
                                  #in zero-matrix

#generating plots

plt.contourf(long,Time,SLP,cmap='twilight_shifted')
plt.xlabel('longitude(degrees)')
plt.ylabel('days since Jan. 1 2015')
plt.title('SLP anomaly (hPa)')
plt.colorbar()
plt.show()

plt.contourf(long,Time,SLP_3,cmap='twilight_shifted')
plt.xlabel('longitude(degrees)')
plt.ylabel('days since Jan. 1 2015')
plt.title('SLP for m = 3 (hPa)')
plt.colorbar()
plt.show()

plt.contourf(long,Time,SLP_5,cmap='twilight_shifted')
plt.xlabel('longitude(degrees)')
plt.ylabel('days since Jan. 1 2015')
plt.title('SLP for m = 5 (hPa)')
plt.colorbar()
plt.show()

