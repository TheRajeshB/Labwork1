"""Created on 2021-09-21

@author: Gabriel Bailey

This code will create a surface plot of u for m=3, n=2, for R<=1, given the z = 11.620.
"""

import numpy as np
from Lab02_IntegrationFunctions import J
import matplotlib.pyplot as plt
from matplotlib import cm
import colorcet as cc
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
fig, ax = plt.subplots(2,2,subplot_kw={"projection": "3d"})

for i in range(2):
    for j in range(2):
        # Plot the surface.
        surf = ax[i][j].plot_surface(X, Y, Z, cmap=cc.cm.bmy,
                        linewidth=0, antialiased=True)
        # Customize the z axis.
        ax[i][j].set_zlim(-1.01, 1.01)
        #ax[i][j].zaxis.set_major_locator(LinearLocator(10))
        # A StrMethodFormatter is used automatically
        ax[i][j].zaxis.set_major_formatter('{x:.02f}')
        ax[i][j].set_xlabel('X')
        ax[i][j].set_ylabel('Y')
        ax[i][j].set_zlabel('Z')


# Add a color bar which maps values to colors.
#fig.colorbar(surf, shrink=0.2, aspect=10)
fig.subplots_adjust(right=0.75)
cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
fig.colorbar(surf, cax=cbar_ax)
fig.suptitle('Surface Plots of Vibrating Membrane')
plt.show()