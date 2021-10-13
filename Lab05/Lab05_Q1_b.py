'''Created on 2021-10-13

@author: Gabriel Bailey

This code will calculate the intesity of light on a screen that has passed through a diffraction grating, and will plot out the results.
'''
import numpy as np
import matplotlib.pyplot as plt
from  cmath  import  pi 

# Define all the constants for this problem (all in m)
width = 20*10**-6
alpha = pi/width
w = 200 * 10**-6
W = 10 * w
lam = 500* 10**-9
f = 1
screen = 10 * 10**-2

N = 5000 # The number of points to try with

# Calculate the intesity using the method provided
xk = lam*f/W * np.arange(N/10)
yk = np.abs(np.sin(alpha*xk))
ck = np.fft.fft(yk, N) #I'm letting numpy pad out the domain for me
I = W**2/N**2 * np.abs(ck)**2
x = lam*f/W * np.arange(N)

# Shift axis to be centered at 0
x = x - lam*f/W*N/2

# Plot the results
fig = plt.figure(figsize=[10,5])
ax = fig.add_subplot(1,1,1)
ax.plot(x, I)
ax.set_xlim( -0.05, 0.05)
ax.set_title('Intensity of Light through Diffraction Grating') 
ax.set_xlabel('x position (m)')
ax.set_ylabel('Intensity')
    
plt.show()
