"""Created on 2021-09-21

@author: Gabriel Bailey

This code will calculate definite integrals using Trapezoidal and Simpsons's rules for integration.
NEED MORE ON WHAT IT PRODUCES
"""
from IntegrationFunctions import trap_integrate, simp_integrate
from time import time
from scipy.constants import pi

import sympy # the symbolic math package

#The function we will be integrating over:
def f(x):
    return(4/(1+x**2))

#ii
print("\npart ii")

print("Real value (via sympy):    ",float(pi))
print("Trapezoidal integration:   ",trap_integrate(f,10,0,1))
print("Simpson's rule integration:",simp_integrate(f,10,0,1))

#iii
print("\npart iii")

#from IntegrationFunctions import symb_integrate
#val = symb_integrate(f,0,1)
val = pi

n = 2
trap_val = trap_integrate(f,n,0,1)
while float(abs(trap_val-val)) > 10**-8: #Multiplicative loop to get to general vicinity fast
    n *= 2
    trap_val = trap_integrate(f,n,0,1)

print("Trap N required to reduce error to O(10^-9):", n)
print("Trap error:", float(abs(trap_val-val)))
start = time()
for i in range(1000):
    trap_integrate(f,n,0,1)
end = time()
print("Trap Time Taken:",(end-start)/1000, "s")

n = 2
simp_val = simp_integrate(f,n,0,1)
while float(abs(simp_val-val)) > 10**-8: #Multiplicative loop to get to general vicinity fast
    #print(float(abs(simp_val-val)))
    n *= 2
    simp_val = simp_integrate(f,n,0,1)

print("Simp N required to reduce error to O(10^-9):", n)
print("Simp error:", float(abs(simp_val-val)))
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
print("The error for N=32 is", e2, "based on the 'practical estimation of errors' method.")
