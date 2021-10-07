# -*- coding: utf-8 -*-
"""
Created on Wed Oct  6 01:06:39 2021

@author: farzh
"""

from Lab04_functions import *
import numpy as np
from cmath import polar
from math import e
import matplotlib.pyplot as plt
#putting in values of constants 

#resistance of resistors in ohms
R1 = 1000 
R3 = 1000
R5 = 1000
R2 = 2000
R4 = 2000
R6 = 2000
 

#capacitance of capacitors in Farads 
C1 =  10**-6 
C2 = 0.5*(10**-6)  

w = 1000 #frequency in rad/s
xplus = 3 #in Volts

#setting up matrix with coefficients of x1, x2 and x3

A = np.array([[(1/R1)+(1/R4)+1j*w*C1 , -1j*w*C1 ,0 ],
             [-1j*w*C1 , (1/R2)+(1/R5)+ 1j*w*C1 + 1j*w*C2, -1j*w*C2],
             [0, -1j*w*C2,(1/R3)+(1/R6)+ 1j*w*C2]],complex)

v = np.array([xplus/R1 , xplus/R2 , xplus/R3],complex) #RHS vector

x = PartialPivot(A,v) #calculating values of x1,x2,x3 using Partial pivot method

print('amp of V1 is :',abs(x[0]),',amp of V2 is:',abs(x[1]),',amp of V3 is:',abs(x[2]),',at t=0')
print('phase of V1 is:',np.angle(x[0]),',phase of V2 is:',np.angle(x[1]),',phase of V3 is:',np.angle(x[2]),',at t=0')

t = np.linspace(0,0.01,100)
#setting V1,V2 and V3 as their respective fuctions of time
V1  = x[0]*(e**(1j*w*t))
V2 = x[1]*(e**(1j*w*t))
V3 = x[2]*(e**(1j*w*t))

plt.plot(t,np.real(V1),label='V1')
plt.plot(t,np.real(V2),label='V2')
plt.plot(t,np.real(V3),label='V3')
plt.xlabel('time in seconds')
plt.ylabel('Voltage (in Volts)')
plt.title('Plots of real voltages vs time for R6 = 2000 Ohm')
plt.legend()
plt.show()


#same code as above , but with R6 = 2000/omega
w = 1000 #frequency in rad/s
R1 = 1000 
R3 = 1000
R5 = 1000
R2 = 2000
R4 = 2000
R6 = 2000/w
 

#capacitance of capacitors in Farads 
C1 =  10**-6 
C2 = 0.5*(10**-6)  


xplus = 3 #in Volts

#setting up matrix with coefficients of x1, x2 and x3

A = np.array([[(1/R1)+(1/R4)+1j*w*C1 , -1j*w*C1 ,0 ],
             [-1j*w*C1 , (1/R2)+(1/R5)+ 1j*w*C1 + 1j*w*C2, -1j*w*C2],
             [0, -1j*w*C2,(1/R3)+(1/R6)+ 1j*w*C2]],complex)

v = np.array([xplus/R1 , xplus/R2 , xplus/R3],complex) #RHS vector

x = PartialPivot(A,v) #calculating values of x1,x2,x3 using Partial pivot method

print('amp of V1 is :',abs(x[0]),',amp of V2 is:',abs(x[1]),',amp of V3 is:',abs(x[2]),',at t=0')
print('phase of V1 is:',np.angle(x[0]),',phase of V2 is:',np.angle(x[1]),',phase of V3 is:',np.angle(x[2]),',at t=0')

t = np.linspace(0,0.01,100)
#setting V1,V2 and V3 as their respective fuctions of time
V1  = x[0]*(e**(1j*w*t))
V2 = x[1]*(e**(1j*w*t))
V3 = x[2]*(e**(1j*w*t))

plt.plot(t,np.real(V1),label='V1')
plt.plot(t,np.real(V2),label='V2')
plt.plot(t,np.real(V3),label='V3')
plt.xlabel('time in seconds')
plt.ylabel('Voltage (in Volts)')
plt.title('Plots of real voltages vs time for R6 = 2 Henry')
plt.legend()
plt.show()