"""Created on 2021-09-21

@author: Gabriel Bailey

This code will calculate definite integrals using Trapezoidal and Simpsons's rules for integration.
"""

import numpy as np
from IntegrationFunctions import J
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator

#Part c
z32 = 11.620 
R=1

# u function for m=3, n=2, t=0
def u32(r,theta,n=2):
    return J(z32*r/R,2)*np.cos(n*theta)

# Make data.
X = np.arange(-R,R,R/100)
Y = np.arange(-R,R,R/100)
X, Y = np.meshgrid(X, Y)

r = np.sqrt(X**2 + Y**2)
theta = np.arctan2(X,Y)
Z = u32(r,theta)
# Mask away outside
Z[r>R] = 0

# Plotting code from https://matplotlib.org/stable/gallery/mplot3d/surface3d.html
fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

# Plot the surface.
surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)

# Customize the z axis.
ax.set_zlim(-1.01, 1.01)
ax.zaxis.set_major_locator(LinearLocator(10))
# A StrMethodFormatter is used automatically
ax.zaxis.set_major_formatter('{x:.02f}')

# Add a color bar which maps values to colors.
fig.colorbar(surf, shrink=0.5, aspect=5)

plt.show()