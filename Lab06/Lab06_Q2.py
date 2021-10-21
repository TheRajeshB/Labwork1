'''Created on 2021-10-20

@author: Gabriel Bailey

This code will simulate the motion of a building using the verlet method plotting the result, and then calculate the normal modes and frequencies.
'''
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from Lab06_functions import adjust_lightness, find_freq_fft, find_freq_lazy

# Generates the matrix for this problem manually, row by row
def generate_A(N):
    A = np.zeros((N,N))
    for i in range(N):
        A[i][i] = -2
        if i > 0:
            A[i][i-1] = 1
        if i < N-1:
            A[i][i+1] = 1
    return A

# Simulates the displacement of floors of a building over time using the Verlet method
def simulate_building(num_floors,maxtime,h,x0,km):
    steps = int(maxtime/h)
    time = np.linspace(0,maxtime,steps) # The time for graphing
    x = np.zeros((2*steps,num_floors)) # The displacement of the floors (includes half steps)
    x[0] = x0
    v = np.zeros((2*steps,num_floors)) # The velocities of the floors (includes half steps)
    A = km*generate_A(num_floors) # The coefficient matrix
    #first half-step
    v[1] = v[0]+1/2*h*A.dot(x[0])

    for t in range(2,2*steps,2):
        x[t] = x[t-2]+h*v[t-1]
        k = h*A.dot(x[t])
        v[t] = v[t-1]+1/2*k
        v[t+1] = v[t-1] + k

    return time, x, v

# Plot the results
def plot_building(num_floors,t,x,freq = None):
    fig = plt.figure(figsize=[10,5])
    ax = fig.add_subplot(1,1,1)
    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    for i in range(num_floors):
        ax.axhline(y=10*i, linestyle=':', color = adjust_lightness(colors[i],0.7))
        ax.plot(t, x[::2,i].T+10*i)
        
    ax.grid(which='both', axis='y')
    ml = MultipleLocator(10)
    ax.yaxis.set_minor_locator(ml)
    ax.set_title('Vibrations of a Building with {} Floors\n (dotted line is 0 for each floor)'.format(num_floors)) 
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('x displacement for floor 0 (cm) ')
    if freq:
        ax.set_title('Vibrations of a Building with {} Floors\n (Normal Mode {})'.format(num_floors,freq[2]))
        ax.text(0.1,num_floors*10-15,'Predicted  Frequency: {:.4f} Hz\nSimulated Frequency: {:.4f} Hz'.format(freq[0],freq[1])) 
        ax.set_xlim([0, 1])
    return ax
    
# Define all the constants for this problem
km = 400 #rad s^-2
dt = 1/1000 #s
maxtime = 1 #s

num_floors = 10 
x0 = np.zeros((1,num_floors))
x0[0][0] = 10 # cm

t10,x10,v10 = simulate_building(num_floors,maxtime,dt,x0,km)
plot_building(num_floors,t10,x10)

num_floors = 3
x0 = np.zeros((1,num_floors))
x0[0][0] = 10 # cm

t3,x3,v3 = simulate_building(num_floors,maxtime,dt,x0,km)
plot_building(num_floors,t3,x3)
plt.show()

#Part b
num_floors = 3
A3 = km*generate_A(num_floors)

# Calculate the eigenvalues and eigenvectors
eigvals, eigvecs = np.linalg.eig(A3)
# np docs: "the column v[:,i] is the eigenvector corresponding to the eigenvalue w[i]"
# Order them in terms of frequency(energy?) (from the eigenvalues)
eigorder = np.argsort(abs(eigvals)) # An array of indices such that eigvals[eigorder] yields a sorted eigvals
eigvals = eigvals[eigorder]
eigvecs = eigvecs[:,eigorder] # Apply it to the column order of eigvecs

# eigenvals = -(angular frequency)^2 = -(2*pi*frequency)^2
freq_predicted = np.sqrt(-eigvals)/(2*np.pi)

# Simulate the 3 normal modes
for i in range(3):
    maxtime = 100
    x0 = eigvecs[:,i]
    t3,x3,v3 = simulate_building(num_floors,maxtime,dt,x0,km)
    freq_lazy = find_freq_lazy(t3,x3[:,0])
    freq_fft = find_freq_fft(maxtime,x3[:,0])
    ax = plot_building(num_floors,t3,x3,[freq_predicted[i],freq_fft, i])
    
    print("Frequency for {}: {}\n Measured {:.4f} or {:.4f}".format(i, freq_predicted[i], freq_lazy, freq_fft))

plt.show()
