# -*- coding: utf-8 -*-
"""
Created on Fri Oct 15 00:22:42 2021

@author: farzh
"""

import numpy as np
import matplotlib.pyplot as plt 

blur_data =np.loadtxt('blur.txt')


def Gaussian(sigma):
    '''Function takes in sigma value and gives point spread function 
    (assuming it is Gaussian distribution) of pixels at positions (x,y)'''
    Matrix = np.zeros(blur_data.shape)#zero matrix with dimensions of original data matrix
    
    for x in range(blur_data.shape[0]): 
        x1=x
        if x1 > (blur_data.shape[0])/2:
            
            x1 -= blur_data.shape[0] #sets rows to right of center to negative values to account for repition
                                     #at end of interval
        for y in range(blur_data.shape[1]):
                y2 = y
                if y2> (blur_data.shape[1])/2:
                  
                  y2 -=  blur_data.shape[1]#sets colums to below of center to negative values to account for repition
                Matrix[x,y]=np.exp(-(x1**2+y2**2)/(2.0*sigma**2)) #computing point spread function
    return Matrix



def unblur():
    '''Function performs Fourier transform of data and point-spread function,
    divides former by latter and takes inverse fourier transform of result which is then divided
    by necessary constant ( depending upon dimensions of matrices involved) to give
    data for unblurred image'''
    sigma = 25 #setting value of sigma to be used in Gaussian/point source function
    Gauss_fourier = np.fft.rfft2(Gaussian(sigma)) #Fourier transform of point-source function
    blur_fourier  = np.fft.rfft2(blur_data)#Fourier transform of image data
    
    blur_Gauss_div = np.zeros(Gauss_fourier.shape,complex) 
    
    for i in range(blur_Gauss_div.shape[0]):
        for j in range(blur_Gauss_div.shape[1]):
            if abs(Gauss_fourier[i,j]) > 10**-3: #check to avoid division by small numbers and zero
                blur_Gauss_div[i,j]= blur_fourier[i,j]/Gauss_fourier[i,j]
            else:
                blur_Gauss_div[i,j]= blur_fourier[i,j]
    Final = np.fft.irfft2(blur_Gauss_div)/((blur_Gauss_div.shape[0]*blur_Gauss_div.shape[1])) #taking inverse Fourier transform
    return Final
#generating plots
plt.imshow(blur_data,cmap='gray')
plt.title('Image from raw data')
plt.show()
plt.imshow(Gaussian(25),vmax=0.001,cmap='gray')
plt.title('Gaussian distribution denisty plot')
plt.show()
plt.imshow(unblur(),cmap='gray')
plt.title('Image after unblurring')
plt.show

