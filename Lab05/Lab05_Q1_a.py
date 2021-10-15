'''Created on 2021-10-13

@author: Gabriel Bailey

This code will calculate the discrete fourier transform of a square wave, sawtooth wave, and modulated sin wave and plot out the resulting frequency domain. It will also plot the functions themselves to confirm that they are the correct ones.
'''
import numpy as np
import matplotlib.pyplot as plt
from  cmath  import  exp,pi 

# Constants
N = 1000
x = np.linspace(0,N,N)

# Generates a 1-cycle square wave
def squarewave(x):
    return np.where(x < N//2,-1,1)

# Generates 1(?) cycle of a sawtooth wave
def sawtooth(x):
    return x

# Generates a modulated sine wave
def mod_sin(x):
    X = len(x)
    return np.sin(pi*x/X)*np.sin(20*pi*x/X)

# Plot the results
def result_plot(x,y,f,c,name):
    fig_test = plt.figure()
    ax_test = fig_test.add_subplot(1,1,1)
    ax_test.plot(x,y)
    ax_test.grid()
    ax_test.set_title(name)
    ax_test.set_xlabel('Time (s)')
    ax_test.set_ylabel('Amplitude')

    bot = np.zeros(len(f))
    fig = plt.figure(figsize=[10,5])
    ax = fig.add_subplot(1,1,1)
    ax.plot(f, c)
    ax.fill_between(f, bot, c)
    ax.set_title('DFT of '+name) 
    ax.set_xlabel('Frequency (Hz)')
    ax.set_ylabel('Amplitude')
     
    plt.show()

#Part a a
y = squarewave(x)
c = abs(np.fft.rfft(y))
f = np.fft.rfftfreq(N)
result_plot(x,y,f,c,"Square Wave")

#Part a b
y = sawtooth(x)
c = abs(np.fft.rfft(y))
f = np.fft.rfftfreq(N)
result_plot(x,y,f,c,"Sawtooth Wave")

#Part a c
y = mod_sin(x)
c = abs(np.fft.rfft(y))
f = np.fft.rfftfreq(N)
result_plot(x,y,f,c,"Modulated Sine Wave")
