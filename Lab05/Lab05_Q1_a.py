'''Created on 2021-10-13

@author: Gabriel Bailey

This code will calculate the discrete fourier transform of a square wave, sawtooth wave, and modulated sin wave and plot out the resulting frequency domain. It will also plot the functions themselves to confirm that they are the correct ones.
'''
import numpy as np
import matplotlib.pyplot as plt
from  numpy import  zeros 
from  cmath  import  exp,pi 

# From textbooks dft.py
def  dft(y): 
    N =  len(y) 
    c  =  zeros(N//2+1,complex) 
    for k  in  range(N//2+1): 
        for  n  in  range(N): 
            c[k] += y[n]*exp(-2j*pi*k*n/N) 
    return  c

# Generates a 1-cycle square wave
def squarewave(x):
    y = np.ones(len(x))
    for i in range(len(x)//2,len(x)):
        y[i] = -1
    return y

# Generates 1(?) cycle of a sawtooth wave
def sawtooth(x):
    return x

# Generates a modulated sine wave
def mod_sin(x):
    X = len(x)
    return np.sin(pi*x/X)*np.sin(20*pi*x/X)

# Plot the results
def result_plot(x,y,c,name):
    fig_test = plt.figure()
    ax_test = fig_test.add_subplot(1,1,1)
    ax_test.plot(x,y)
    ax_test.grid()
    ax_test.set_title(name)
    ax_test.set_xlabel('Time (s)')
    ax_test.set_ylabel('Amplitude')
    
    fig = plt.figure(figsize=[10,5])
    ax = fig.add_subplot(1,1,1)
    ax.bar(range(len(c)), c, width=1.3)
    ax.set_title('DFT of '+name) 
    ax.set_xlabel('Frequency (Hz)')
    ax.set_ylabel('Amplitude')
     
    plt.show()

#Part a a
N = 1000
x = np.linspace(0,N,N)

y = squarewave(x)
c = abs(dft(y))
result_plot(x,y,c,"Sqaure Wave")

#Part a b
y = sawtooth(x)
c = abs(dft(y))
result_plot(x,y,c,"Sawtooth Wave")

#Part a c
y = mod_sin(x)
c = abs(dft(y))
result_plot(x,y,c,"Modulated Sine Wave")
