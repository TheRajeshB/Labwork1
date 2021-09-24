from math import *
import numpy as np
from scipy.integrate import quad
from scipy.special import kn
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator
      

def V(r,z,u):
    Q= 10**-13
    m = tan(u)
    n = cos(u)
    epsilon = 8.854e-12 #epsilon-nought in calculation of electrostatic constant
    l = 0.001 #l in mm
    k = Q/(4*pi*epsilon)
    c = ((z-l*m)**2 + r**2)**0.5 #terms under root in denominator
        
    return (k*e**-(m**2))/((n**2)*c)
    
    
#fix to simpsons
 
def simp(f,a,b,N):
    

     h= (b-a)/N
    
     sum_odd = 0.0
     sum_even = 0.0
     initial_sum = f(a)+f(b)
    
     trap_area=[]
    
    #calculating last slice in trapezoidal if N is odd 
     if N % 2 == 1: #if N is odd calc the last slice with trapezoidal rule
         trap_area= np.append(trap_area,h*f(a+(N-1)*h))
         N -= 1
         b= b-h
         #loop to find sum of odd terms for Simpson integration
     for k in range(1,int(N),2):
        sum_odd  += f((a+k*h))
     
    #loop to find sum of even terms for Simpson integration
     for k in range(2,int(N),2):
        sum_even += f((a+k*h))
     
     return(h*(initial_sum + 4*sum_odd + 2*sum_even)/3) + (trap_area)
     
     


#print((simp(V,-pi/2,pi/2,9)))
from IntegrationFunctions import trap_integrate, simp_integrate

R = 0.005
# Make data.
z = np.linspace(-0.005,0.005,100)#range of values for z in 'm'
r = np.linspace(0.00025,0.005,100)#range of values for r
r, z = np.meshgrid(r, z)

v = simp_integrate(lambda u : V(r,z,u),9,-pi/2,pi/2,)

# Plotting code from https://matplotlib.org/stable/gallery/mplot3d/surface3d.html
fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

# Plot the surface.
surf = ax.plot_surface(r, z, v, cmap=cm.coolwarm,
                       linewidth=1, antialiased=False)

# Customize the z axis.
ax.set_zlim(-1.01, 1.01)
ax.zaxis.set_major_locator(LinearLocator(10))
# A StrMethodFormatter is used automatically
ax.zaxis.set_major_formatter('{x:.02f}')
ax.autoscale()

# Add a color bar which maps values to colors.
fig.colorbar(surf, shrink=0.5, aspect=5)

plt.show()