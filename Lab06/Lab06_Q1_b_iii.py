# -*- coding: utf-8 -*-
"""
Created on Thu Oct 21 01:24:16 2021

@author: farzh
"""

import numpy as np 
import matplotlib.pyplot as plt

#defining constants
w_0 = 1 
v_p = 1
tau = 1 #tau
g = 0.5 #gamma
v_f = 0.1
#parameters for RK4 method
N = 5000
t_start = 0 
t_end = 300
h = (t_end-t_start)/N

#initial conditions
x_initial = 1.000227#intial condition for position
y_initial = -1#initial velocity 
def f(r,t):
    x = r[0]
    y = r[1]  #defining y = dx/dt
    f_x = y #ODE for velocity 
    f_y = -(w_0*(x - v_p*t) + (abs(y)/tau) + g*np.exp(-abs(y)/v_f)) #defining ODE for equation of motion
    return np.array([f_x,f_y],float)

tpoints = np.arange(t_start,t_end,h) #points to be considered of each time-step
r = np.array([x_initial,y_initial],float)
x_values=[] #list to store values of x(t)
y_values=[] #list to store values of dx/dt

#solving ODE using RK4 method
for t in tpoints:
  x_values.append(r[0])  
  y_values.append(r[1])
  k1 = h*f(r, t)
  k2 = h*f(r+0.5*k1, t+0.5*h)
  k3 = h*f(r+0.5*k2, t+0.5*h)
  k4 = h*f(r+k3, t+h)
  r += (k1 + 2*k2 + 2*k3 + k4)/6
 
x_array = np.array([x_values])  
y_array = np.array([y_values])  
#energy = ((w_0**2)*(x_array - v_p*tpoints)**2) + y_array**2

#generating plots

#plt.plot(tpoints,energy[0])
#plt.title('Graph showing decay of energy')
#plt.ylabel('energy (Joule)')
#plt.xlabel('time (second)')
#plt.show()
  
plt.plot(tpoints,x_values)
plt.xlabel('time (second)')
plt.ylabel('x(t) (meter)')
plt.title('plot for solution of equation of motion ')
plt.show()

plt.plot(tpoints,y_values)
plt.xlabel('time (seconds)')
plt.ylabel('dx/dt (m/s)')
plt.title('graph for dx/dt vs time')
plt.show()


#same code as above but with values adjusted to solve equation from Q1(a)v.
w_0 = 1 
v_p = 0.05
tau = 1 #tau
g = 0.5 #gamma
v_f = 0.1
#parameters for RK4 method
N = 5000
t_start = 0 
t_end = 300
h = (t_end-t_start)/N

x_initial = 1.0000#intial condition
y_initial = 1
def f(r,t):
    x = r[0]
    y = r[1]  #defining y = dx/dt
    f_x = y #ODE for velocity 
    f_y = -(w_0*(x - v_p*t) + (abs(y)/tau) + g*np.exp(-abs(y)/v_f)) #defining ODE for equation of motion
    return np.array([f_x,f_y],float)

tpoints = np.arange(t_start,t_end,h) #points to be considered of each time-step
r = np.array([x_initial,y_initial],float)
x_values=[] #list to store values of x(t)
y_values=[] #list to store values of dx/dt


for t in tpoints:
  x_values.append(r[0])  
  y_values.append(r[1])
  k1 = h*f(r, t)
  k2 = h*f(r+0.5*k1, t+0.5*h)
  k3 = h*f(r+0.5*k2, t+0.5*h)
  k4 = h*f(r+k3, t+h)
  r += (k1 + 2*k2 + 2*k3 + k4)/6

plt.plot(tpoints,x_values)
plt.xlabel('time (second)')
plt.ylabel('x(t) (meter)')
plt.title('plot for solution of equation of motion ')
plt.show()

plt.plot(tpoints,y_values)
plt.xlabel('time (seconds)')
plt.ylabel('dx/dt (m/s)')
plt.title('graph for dx/dt vs time')
plt.show()
