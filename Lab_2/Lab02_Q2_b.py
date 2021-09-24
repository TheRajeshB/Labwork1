"""Created on 2021-09-21

@author: Gabriel Bailey

This code will calculate the 0th, 3rd, and 5th order Bessel functions between 0 and 20 using Simpson's rule with N=1000. It will then plot the results, and the delta of the results to the scipy versions of the functions.
"""

import numpy as np
from scipy.special import jv, kn
from IntegrationFunctions import J # Where the Bessel function is defined
import matplotlib.pyplot as plt

#Create the x values
xs = np.linspace(0, 20, 1000)

# Calculate the y values using my Bessel function
j0 = np.array([J(x,0) for x in xs])
j3 = np.array([J(x,3) for x in xs])
j5 = np.array([J(x,5) for x in xs])

# Calculate the y values using scipy's Bessel function
j0sp = np.array([jv(0,x) for x in xs])
j3sp = np.array([jv(3,x) for x in xs])
j5sp = np.array([jv(5,x) for x in xs])


#Plot my Bessel function values
plt.figure()
plt.errorbar(xs, j0, label = "J0")
plt.errorbar(xs, j3, label = "J3")
plt.errorbar(xs, j5, label = "J5")
plt.legend()
plt.xlabel('x')
plt.ylabel('Bessel Function Value')
plt.title('Plot of Our Bessel Functions')  

#Plot absolute difference from scipy Bessel function values
plt.figure()
plt.errorbar(xs, abs(j0sp-j0), label = "J0")
plt.errorbar(xs, abs(j3sp-j3), label = "J3")
plt.errorbar(xs, abs(j5sp-j5), label = "J5")
plt.legend()
plt.xlabel('x')
plt.ylabel('Bessel Function Value')
plt.title('Absolute Difference from\n Scipy Bessel Functions')  

plt.show()