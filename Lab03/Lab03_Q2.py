"""Created on 2021-09-29

@author: Gabriel Bailey

This code will calculate the period of a spring using gaussian quadrature with several values of N, and will then print out the results and there errors as well as create 3 plots: one of what each sample point evaluates to, one of those values after weighting, and on that compares the results of our calculation with classical and relativistic limits.
"""

import numpy as np
#import sympy as sp
from scipy.constants import c, pi
from functions_Lab03 import gaussxw, gaussxwab, gaussxwab_convert
import matplotlib.pyplot as plt

def g(x,x0,k,m):
    #print(k*(x0**2-x**2)*(2*m*c**2 + k*(x0**2-x**2)/2), 2*m*c**2)
    return c*( k*(x0**2-x**2)*(2*m*c**2 + k*(x0**2-x**2)/2) / (2*(m*c**2 + k*(x0**2-x**2)/2)**2) )**(1/2)

def T_details(x0, k, m, N):
    x,w = gaussxwab(N,0,x0)
    gk, wgk = np.empty(N), np.empty(N)
    sum = 0.0
    for i in range(N):
        gk[i] = 4/g(x[i],x0,k,m)
        wgk[i] = w[i]*gk[i]
        sum += wgk[i]
    return sum, x, gk, wgk

x0 = 0.01   # m
m = 1       # kg
k = 12      # N/m
T8, x8, gk8, wgk8 = T_details(x0,k,m,8)
T16, x16, gk16, wgk16 = T_details(x0,k,m,16)
T32 = T_details(x0,k,m,32)[0]
T_class = 2*pi*(m/k)**(1/2)


# Testing our Gaussian quadrature method against our other methods:

#from functions_Lab03 import symb_integrate, trap_integrate, simp_integrate
# def g_int(x):
#     return 4/g(x,x0,k,m)
# x00 = 10**(-15)

# Ttrap = trap_integrate(g_int,100000000, 0, x0-x00)
# print(Ttrap)
# 1.8265933061561042 for N = 100000000

# Tsimp = simp_integrate(g_int,100000000, 0, x0-x00)
# print(Tsimp)
# 1.8223053758301062 for N = 100000000

# Tsymb = symb_integrate(g_int, 0, x0) doesn't work at all
# print(Tsymb)

#2.0419674852319942**-17 value from wolfram alpha
#-0.023094 from wa 2nd try
#0.023430709023252595 from sympy

# Fractional error estimation
fra_err8 = abs((T16-T8)/T8)
fra_err16 = abs((T32-T16)/T16)

# Print the results
print("Period (classical):", T_class,"s")
print("Period with N = 8 :",T8,"s")
print("Relative error estimated   :",fra_err8)
print("Relative error vs classical:",abs((T8-T_class)/T_class))
print("Period with N = 16:",T16,"s")
print("Relative error estimated   :",fra_err16)
print("Relative error vs classical:",abs((T16-T_class)/T_class))


# Plot the results for N=8 and N=16
fig1 = plt.figure(figsize=[10,5])
ax = fig1.add_subplot(1,1,1)
ax.plot(x8, gk8, 'o-', label = "N = 8")
ax.plot(x16, gk16, 's-', label = "N = 16")
ax.legend()
ax.grid()
ax.set_xlabel('x (m)')
ax.set_ylabel('Value of 4/g_k (s/m)')
ax.set_title('Plot of 4/g_k for x0=1cm, m=1kg, k=12N/m')  

fig2 = plt.figure(figsize=[10,5])
ax = fig2.add_subplot(1,1,1)
ax.plot(x8, wgk8, 'o-', label = "N = 8")
ax.plot(x16, wgk16, 's-', label = "N = 16")
ax.legend()
ax.grid()
ax.set_xlabel('x (m)')
ax.set_ylabel('Value of 4w_k/g_k (s/m)')
ax.set_title('Plot of Weighted Values for 4w_k/g_k for x0=1cm, m=1kg, k=12N/m')  

plt.show()

# Start of part b

xc=m*k*c**2

# Start of part c

T200, x200, gk200, wgk200 = T_details(x0,k,m,200)
T400 = T_details(x0,k,m,400)[0]
fra_err200 = abs((T400-T200)/T200)

print("Period with N = 200:",T200,"s")
print("Relative error estimated   :",fra_err200)

N = 200
x_gen, weights = gaussxw(N)

def T(x0):
    x, w = gaussxwab_convert(x_gen,weights,0,x0)
    sum = 0.0
    for i in range(N):
        sum += w[i]*4/g(x[i],x0,k,m)
    return sum

m = 10 # Number of points on plot
xs = np.linspace(1,xc,m)
Ts = np.empty(m)
for i in range(m):
    Ts[i] = T(xs[i])

fig3 = plt.figure(figsize=[10,5])
ax = fig3.add_subplot(1,1,1)
ax.plot(xs, Ts, 'o-', label = "Calculated Period")
plt.axhline(y=T_class, color='g', linestyle='-', label = "Classical Limit") # Classical limit
ax.plot([1,xc], [4*1*c,4*xc/c], '-', label = "Relativistic Limit")
ax.legend()
ax.grid()
ax.set_xlabel('x0 (m)')
ax.set_ylabel('Period T (s)')
ax.set_title('Period with Different Initial Values for x')  

plt.show()