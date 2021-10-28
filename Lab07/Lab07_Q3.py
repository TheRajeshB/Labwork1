# Adapted from bulirsch.py
"""Created on 2021-10-27

@author: Gabriel Bailey

This code adapts the code from bulirsch.py to work on an adaptive timestep, and simulates the Belousov--Zlzabotinsky  reaction, creating a plot of the chemical concentrations.
"""
import numpy as np
import matplotlib.pyplot as plt

# Constants
a = 1
b = 3
x0 = 0
y0 = 0
delta = 1e-10     # Required position accuracy per unit time
n_max = 8

# recursion step
# else recall r t h, break it in half pg 3
def f(con):
    x = con[0]
    y = con[1]
    vx = 1-(b+1)*x + a*x**2*y
    vy = b*x - a*x**2*y
    return np.array([vx,vy],float)

def step(r, t, H):
    # The actual arrays we are trying to fill with data
    tpoints = []
    rpoints = []

    # Do one modified midpoint step to get things started
    n = 1
    r1 = r + 0.5*H*f(r)
    r2 = r + H*f(r1)

    R1 = np.empty([1,2],float)
    R1[0] = 0.5*(r1 + r2 + 0.5*H*f(r2))

    # Now increase n until the required accuracy is reached
    error = 2*H*delta
    while error>H*delta:
        n += 1
        if n > n_max:
            # Go to a smaller time step
            t1,r1 = step(r, t, H/2)
            t2,r2 = step(r1[-1], t+H/2, H/2)
            # Place new data into existing lists
            tpoints.extend(t1)
            tpoints.extend(t2)
            rpoints.extend(r1)
            rpoints.extend(r2)
            return tpoints,rpoints
        h = H/n

        # Modified midpoint method
        r1 = r + 0.5*h*f(r)
        r2 = r + h*f(r1)
        for i in range(n-1):
            r1 += h*f(r2)
            r2 += h*f(r1)

        # Calculate extrapolation estimates.  Arrays R1 and R2
        # hold the two most recent lines of the table
        R2 = R1
        R1 = np.empty([n,2],float)

        R1[0] = 0.5*(r1 + r2 + 0.5*h*f(r2))
        for m in range(1,n):
            epsilon = (R1[m-1]-R2[m-1])/((n/(n-1))**(2*m)-1)
            R1[m] = R1[m-1] + epsilon
        if np.isnan(epsilon[0]) or np.isnan(epsilon[1]):
            epsilon = np.array([np.inf,np.inf])
        error = abs(epsilon[0])

    # Set r equal to the most accurate estimate we have,
    # before moving on to the next big step
    r = R1[n-1]
    # Add current point to data lists
    tpoints.append(t+H)
    rpoints.append(r)
    return tpoints,rpoints

# Initial simulation conditions
r0 = np.array([x0,y0],float)
t0 = 0.0
H = 20 # Size of inital step, also total size

# Simulate
tp,rp = step(r0, t0, H)

# Format for plotting
tpoints = [t0] + tp
concentration= [r0] + rp
tpoints, concentration = np.array(tpoints), np.array(concentration)

# size of h
dtpoints = np.array(tpoints[1:]) - np.array(tpoints[:-1])

# Plot the results
plt.figure()
plt.plot(tpoints,concentration[:,0],".-", label="x")
plt.plot(tpoints,concentration[:,1],".-", label="y")
plt.plot(tpoints[:-1], dtpoints, label = "Timestep Size (s)") # drop the last point in tpoints
plt.xlabel("Time (s)")
plt.ylabel("Chemical concentrations")
plt.title('Belousov-Zhabotinsky Reaction over Time')
plt.legend()
plt.grid()
plt.tight_layout()
plt.savefig('Q3.png', dpi=150)

plt.show()
