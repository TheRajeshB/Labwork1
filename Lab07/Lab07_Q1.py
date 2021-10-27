# Solution to Newman 8.8, Space garbage.
# Author: Nico Grisouard, Univ. of Toronto
"""Created on 2021-10-27

@author: Gabriel Bailey

This code adapts the code from above to repeat the calculation using adaptive timestepping, and outputs a comparison plot, prints the time they take, and a plot of the adaptive step size over time.
"""


import numpy as np
import matplotlib.pyplot as plt
from time import time


def rhs(r):
    """ The right-hand-side of the equations
    INPUT:
    r = [x, vx, y, vy] are floats (not arrays)
    note: no explicit dependence on time
    OUTPUT:
    1x2 numpy array, rhs[0] is for x, rhs[1] is for vx, etc"""
    M = 10.
    L = 2.

    x = r[0]
    vx = r[1]
    y = r[2]
    vy = r[3]

    r2 = x**2 + y**2
    Fx, Fy = - M * np.array([x, y], float) / (r2 * np.sqrt(r2 + .25*L**2))
    return np.array([vx, Fx, vy, Fy], float)


# %% This next part adapted from Newman's odesim.py --------------------------|
a = 0.0
b = 10.0
N = 10000  # let's leave it at that for now
h = (b-a)/N
start_na = time()
tpoints = np.arange(a, b, h)
xpoints = []
vxpoints = []  # the future dx/dt
ypoints = []
vypoints = []  # the future dy/dt

# below: ordering is x, dx/dt, y, dy/dt
r = np.array([1., 0., 0., 1.], float)
for t in tpoints:
    xpoints.append(r[0])
    vxpoints.append(r[1])
    ypoints.append(r[2])
    vypoints.append(r[3])
    k1 = h*rhs(r)  # all the k's are vectors
    k2 = h*rhs(r + 0.5*k1)  # note: no explicit dependence on time of the RHSs
    k3 = h*rhs(r + 0.5*k2)
    k4 = h*rhs(r + k3)
    r += (k1 + 2*k2 + 2*k3 + k4)/6
end_na = time()


h = 0.01
delta = 10**-6

start_ad = time()
# below: ordering is x, dx/dt, y, dy/dt
r = np.array([1., 0., 0., 1.], float)

t_ad = [0]
x_ad = [r[0]]
vx_ad = [r[1]]  # the future dx/dt
y_ad = [r[2]]
vy_ad = [r[3]]  # the future dy/dt

#
def take_step(r,h):
    
    k1 = h*rhs(r)  # all the k's are vectors
    k2 = h*rhs(r + 0.5*k1)  # note: no explicit dependence on time of the RHSs
    k3 = h*rhs(r + 0.5*k2)
    k4 = h*rhs(r + k3)
    return r + (k1 + 2*k2 + 2*k3 + k4)/6

while t_ad[-1] < b:
    h = min(h, (b-t_ad[-1])/2) # Make it end at b
    r0 = take_step(r,h)
    r1 = take_step(r0,h)
    r2 = take_step(r, 2*h)
    
    rho = h*delta/(np.sqrt((r1[0]-r2[0])**2+(r1[2]-r2[2])**2)/30)
    
    if rho >= 1:
        r = r1
        t_ad.append(t_ad[-1]+2*h)
        h = h*max(rho**(1/4),2) # Prevent diverging result
    else:
        h = h*rho**(1/4)
        r0 = take_step(r,h)
        r1 = take_step(r0,h)
        r = r1
        t_ad.append(t_ad[-1]+2*h)

    x_ad.append(r[0])
    vx_ad.append(r[1])
    y_ad.append(r[2])
    vy_ad.append(r[3])
end_ad = time()
print("Time taken for nonadaptive timestep method: {} s".format(end_na-start_na))
print("Time taken for    adaptive timestep method: {} s".format(end_ad-start_ad))
#print(x_ad, y_ad)
plt.figure()
plt.plot(xpoints, ypoints, ':', label = "Nonadaptive")
plt.plot(x_ad, y_ad, '.', label = "Adaptive")
plt.xlabel("$x$")
plt.ylabel("$y$")
plt.legend()
plt.title('Trajectory of a ball bearing around a space rod.')
plt.axis('equal')
plt.grid()
plt.tight_layout()
plt.savefig('Q1a.png', dpi=150)

plt.figure()
dt_ad = np.array(t_ad[1:]) - np.array(t_ad[:-1])
plt.plot(t_ad[:-1], dt_ad, label = "Timestep") # drop the last point in tpoints
r_ad =  np.dstack((np.array(x_ad),np.array(y_ad)))
r_ad = np.array([np.linalg.norm(x) for x in r_ad[0]])
r_ad *= np.max(dt_ad)
plt.plot(t_ad[:-1], r_ad[:-1], c="lightblue" , label="Distance from rod (scaled)")
plt.xlabel("Time t (s)")
plt.ylabel("Timestep h (s)")
plt.title('Size of Timestep over Time')
plt.legend()
plt.grid()
plt.tight_layout()
plt.savefig('Q1c.png', dpi=150)

plt.show()