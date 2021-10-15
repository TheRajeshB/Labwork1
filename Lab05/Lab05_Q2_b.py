"""Created on 2021-10-14

@author: Gabriel Bailey

This code will process an audio file to generate a plot of the sound, a plot of the sound after filtering, and save the modified file.
"""

import numpy as np
from scipy.constants import c, pi
import matplotlib.pyplot as plt
from scipy.io.wavfile import read, write
from numpy import empty
import bisect

def create_time_plot(t,channel_0,channel_1,filtered =''):
    fig, axs = plt.subplots(2, 1, constrained_layout=True)
    for i in range(2):
        axs[i].plot(t, [channel_0,channel_1][i])
        axs[i].set_title('Channel '+str(i))
        axs[i].set_xlabel('Time (s)')
        axs[i].set_ylabel('Sound Pressure')

    fig.suptitle('Air Pressure of GraviteaTime.wav over Time '+filtered, fontsize=16)

def create_freq_plot(f,channel_0,channel_1,filtered =''):
    fig, axs = plt.subplots(2, 1, constrained_layout=True)
    for i in range(2):
        axs[i].plot(f, [channel_0,channel_1][i])
        axs[i].set_title('Channel '+str(i))
        axs[i].set_xlabel('Frequency (Hz)')
        axs[i].set_ylabel('Amplitude')

    fig.suptitle('Frequency Distribution of GraviteaTime.wav '+filtered, fontsize=16)


# read the data
# sample is the sampling rate, data is the data in each channel
sample, data = read('GraviteaTime.wav')
# sample is the sampling frequency, 44100 Hz
# separate into channels
channel_0 = data[:, 0]
channel_1 = data[:, 1]
N = len(channel_0)

timestep = 1/sample
time = N/sample

t = np.linspace(0,time,N)
create_time_plot(t,channel_0,channel_1)

f = np.fft.rfftfreq(N,timestep)
c0 = np.fft.rfft(channel_0)
c1 = np.fft.rfft(channel_1)
#create_freq_plot(f,c0,c1)

# Just lazily finding the index of 880Hz with binary search
cutoff = bisect.bisect_left(f, 880)
c0[cutoff:] = 0
c1[cutoff:] = 0

#create_freq_plot(f,c0,c1)
channel_0_out = np.fft.irfft(c0)
channel_1_out = np.fft.irfft(c1)

create_time_plot(t,channel_0_out,channel_1_out,'Filtered')
plt.show()

# this creates an empty array data_out with the same shape as "data"
# (2 x N_Points) and the same type as "data" (int16)
data_out = empty(data.shape, dtype = data.dtype)
# fill data_out
data_out[:, 0] = channel_0_out
data_out[:, 1] = channel_1_out
write('GraviteaTime_lpf.wav', sample, data_out)

print("Done")

