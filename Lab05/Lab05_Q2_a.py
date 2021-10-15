"""Created on 2021-10-13

@author: Gabriel Bailey

This code will calculate the period of a spring using gaussian quadrature with several values of N, and will then print out the results and there errors as well as create 3 plots: one of what each sample point evaluates to, one of those values after weighting, and on that compares the results of our calculation with classical and relativistic limits.
"""

import numpy as np
from scipy.constants import c, pi
import matplotlib.pyplot as plt
from Lab05_functions import T_details

# Set constants
m = 1       # kg
k = 12      # N/m
xc=np.sqrt(m/k)*c
N = [5000,100000,300000]
time = [50,50,150]
x0 = [1,xc,10*xc]
labels = ["1 m", "x_c", "10x_c"]
colors = [['blue','darkblue'],['orange','darkorange'],['green','darkgreen']]

# Calculates position and time for a particle on a spring with starting position x0 at N points from 0 to 'time' seconds.
def euler_cromer(x0,N,time):
    x = np.empty(N)
    v = np.empty(N)
    x[0] = x0
    v[0] = 0
    delta_t = time/N
    for i in range(1,N):
        v[i] = v[i-1] - k/m*x[i-1]*np.power(1-v[i-1]**2/c**2, 3/2) * delta_t
        x[i] = x[i-1] + v[i] * delta_t
        #if x[i-1] > -100:
        #    print(i,x[i-1],x[i],v[i-1],v[i],np.power(1-v[i-1]**2/c**2, 3/2))
    return x, v

# Adds a time-domain plot. Last 3 arguments are titles/formatting
def add_time_plot(t, x, v, N, time, startname, case):
    fig_time = plt.figure(figsize=[10,6])
    ax_time = fig_time.add_subplot(1,1,1)
    ax_time.plot(t,x, label = 'Position (m)')
    ax_time.plot(t,v, label = 'Velocity (m/s)')
    ax_time.set_xlabel('Time (s)')
    ax_time.set_ylabel('Position (m)\nVelocity (m/s)')
    #ax_vis.set_yscale('log')
    ax_time.set_title('Position of Particle on Relativistic\n Spring over Time (Case {}:\nx_0 = {}, N = {}, dt = {}'.format(case,startname,N,time/N))
    ax_time.legend()#title="Starting position")

# Adds a frequency-domain plot. Last 4 arguments are titles/formatting
def add_freq_plot(ax_freq, f, c, T_expected, N, time, startname, case, colors, title):
    maxc = np.max(c)
    ax_freq.plot(f[:N//2],abs(c[:N//2]/maxc), label = 'Case {}: x_0 = {}, N = {}, dt = {}'.format(case,startname,N,time/N),c=colors[0])
    ax_freq.axvline(x=1/T_expected, label='Case {} Predicted Frequency'.format(case), c=colors[1],ls = '--')
    ax_freq.set_xlabel('Frequency (Hz)')
    ax_freq.set_xlim(0,5)
    ax_freq.set_ylabel('Magnitude (normalized to max of 1)')
    #ax_vis.set_yscale('log')
    ax_freq.set_title('Frequency of '+title+' of Particle on Relativistic Spring')
    ax_freq.legend()#title="Starting position")
    
# Prints the relative difference due to the non-stability of the system over the course of the simulation
def print_errors(x,v,n,i):
    print("Case {} nonstability:".format(i))
    print("Position relative error: ",abs(np.max(x[-n//10:])-np.max(x[:n//10]))/np.max(x[:n//10]))
    print("Velocity relative error: ",abs(np.max(v[-n//10:])-np.max(v[:n//10]))/np.max(v[:n//10]))

# Initialize nested lists to store the 3 cases:
t =[[0],[0],[0]]
x= [[0],[0],[0]]
v= [[0],[0],[0]]
cx= [[0],[0],[0]]
cv= [[0],[0],[0]]
f= [[0],[0],[0]]
T_predicted =[[0],[0],[0]]

for i in range(3):
    t[i] = np.linspace(0,time[i],N[i])
    x[i],v[i] = euler_cromer(x0[i], N[i], time[i])
    # Uncomment to see time plots:
    # print_errors(x[i],v[i],n[i],i+1)
    # add_time_plot(t[i], x[i], v[i], N[i], time[i], labels[i], i+1)
# plt.show()

fig_freqx = plt.figure(figsize=[10,6])
ax_freqx = fig_freqx.add_subplot(1,1,1)
for i in range(0,3):
    cx[i] = np.fft.rfft(x[i])
    f[i] = np.fft.fftfreq(N[i],time[i]/N[i])
    T_predicted[i] = T_details(x0[i],k,m,32)[0]
    add_freq_plot(ax_freqx,f[i],cx[i],T_predicted[i],N[i],time[i],labels[i],i+1,colors[i], "Position")

fig_freqv = plt.figure(figsize=[10,6])
ax_freqv = fig_freqv.add_subplot(1,1,1)
for i in range(0,3):
    cv[i] = np.fft.rfft(v[i])
    f[i] = np.fft.fftfreq(N[i],time[i]/N[i])
    #T_predicted[i] = T_details(x0[i],k,m,32)[0] #Uncomment if above loop is commented out
    add_freq_plot(ax_freqv,f[i],cv[i],T_predicted[i],N[i],time[i],labels[i],i+1,colors[i], "Velocity")

plt.show()
print("Done")

