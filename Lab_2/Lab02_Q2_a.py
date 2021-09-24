"""Created on 2021-09-21

@author: Gabriel Bailey

This code will calculate definite integrals using Trapezoidal and Simpsons's rules for integration, and print the results for Q2a part ii to iv in the terminal.
"""

from time import time
from scipy.constants import pi
#import our functions for integration with Trapezoidal and Simpsons's rules:
from IntegrationFunctions import trap_integrate, simp_integrate

#The function we will be integrating over:
def f(x):
    return(4/(1+x**2))

#ii
print("\npart ii")

print("Real value:                ",float(pi))
print("Trapezoidal integration:   ",trap_integrate(f,4,0,1))
print("Simpson's rule integration:",simp_integrate(f,4,0,1))

#iii
print("\npart iii (TA said we could have this as part of the Q2a code)")

# The actual value of the integration
val = pi

# Set initial n and value of trapezoidal integration:
n = 2
trap_val = trap_integrate(f,n,0,1)

# Loop until the error is less than O(10^-8), multiplying n by 2 each time:
while float(abs(trap_val-val)) > 10**-8:
    n *= 2
    trap_val = trap_integrate(f,n,0,1)

#Print results
print("\nTrap N required to reduce error to O(10^-9):", n)
print("Trap error:", float(abs(trap_val-val)))

# Time the function by running it 1000 times and dividing the total run time by 1000
start = time()
for i in range(1000):
    trap_integrate(f,n,0,1)
end = time()
print("Trap Time Taken:",(end-start)/1000, "s")

# Set initial n and value of Simpsons's integration:
n = 2
simp_val = simp_integrate(f,n,0,1)

# Loop until the error is less than O(10^-8), multiplying n by 2 each time:
while float(abs(simp_val-val)) > 10**-8:
    n *= 2
    simp_val = simp_integrate(f,n,0,1)

#Print results
print("\nSimp N required to reduce error to O(10^-9):", n)
print("Simp error:", float(abs(simp_val-val)))

# Time the function by running it 1000 times and dividing the total run time by 1000
start = time()
for i in range(1000):
    trap_integrate(f,n,0,1)
end = time()
print("Simp Time Taken:",(end-start)/1000, "s")

#iv
print("\npart iv")

I1 = trap_integrate(f,16,0,1)
I2 = trap_integrate(f,32,0,1)
e2 = 1/3*abs(I2-I1)
print("The error for N=32 using Simpsons's method:")
print("The 'practical estimation of errors' error is", e2)
print("The actual error is                          ", float(abs(I2-val)))
