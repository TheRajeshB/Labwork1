"""Created on 2021-09-21

@author: Gabriel Bailey

This code will calculate definite integrals using Trapezoidal and Simpsons's rules for integration.
"""

import numpy as np
from scipy.special import jv, kn
from IntegrationFunctions import J
import matplotlib.pyplot as plt

xs = np.arange(0, 20, 20/1000) # Use linspace

j0 = np.array([J(x,0) for x in xs])
j3 = np.array([J(x,3) for x in xs])
j5 = np.array([J(x,5) for x in xs])

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

#Plot scipy Bessel function values
plt.figure()
plt.errorbar(xs, j0sp-j0, label = "J0")
plt.errorbar(xs, j3sp-j3, label = "J3")
plt.errorbar(xs, j5sp-j5, label = "J5")
plt.legend()
plt.xlabel('x')
plt.ylabel('Bessel Function Value')
plt.title('Difference from Scipy Bessel Functions')  

plt.show()