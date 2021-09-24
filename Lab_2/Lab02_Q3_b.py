"""Created on 2021-09-21

@author: Farz Halim
@author: Gabriel Bailey

This code will calculate the electrostatic potential V of an infinite line of charge as a function of radius and z using Simpson's rule. It will then plot it out as a surface plot.
"""

from math import *
import numpy as np
from scipy.integrate import quad
from scipy.special import kn
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator
from IntegrationFunctions import simp_integrate

# The part of V that is inside the integral:
def V(r,z,u):
    Q= 10**-13
    m = tan(u)
    n = cos(u)
    epsilon = 8.854e-12 #epsilon-nought in calculation of electrostatic constant
    l = 0.001 #l in mm
    k = Q/(4*pi*epsilon)
    c = ((z-l*m)**2 + r**2)**0.5 #terms under root in denominator
        
    return (k*e**-(m**2))/((n**2)*c)

R = 0.005
# Make data.
z = np.linspace(-0.005,0.005,100)#range of values for z in 'm'
r = np.linspace(0.00025,0.005,100)#range of values for r
r, z = np.meshgrid(r, z)

v = simp_integrate(lambda u : V(r,z,u),9,-pi/2,pi/2,)

# Plotting code from https://matplotlib.org/stable/gallery/mplot3d/surface3d.html
fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

# Plot the surface.
surf = ax.plot_surface(r, z, v, cmap=cm.rainbow,
                       linewidth=1, antialiased=True)

# Customize the z axis.
ax.set_zlim(-1.01, 1.01)
ax.zaxis.set_major_locator(LinearLocator(10))
# A StrMethodFormatter is used automatically
ax.zaxis.set_major_formatter('{x:.02f}')
ax.autoscale()

ax.set_xlabel('r')
ax.set_ylabel('z')
ax.set_zlabel('V')

# Add a color bar which maps values to colors.
fig.subplots_adjust(right=0.82)
cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
fig.colorbar(surf, cax=cbar_ax)
#fig.colorbar(surf, shrink=0.5, aspect=5)
fig.suptitle('Electrostatic Potential over r and z')
plt.show()