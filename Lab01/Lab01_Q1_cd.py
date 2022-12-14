"""Created on Mon Sep 13 14:06:42 2021

@author: farzh
"""
#1c 

import matplotlib.pyplot as plt
import numpy as np

def Mercury_orbit():
    
    #setting time-step and number of steps for integration
    delta_t = 0.0001
    years = 1.0 #number of Earth years over which values are calculated
    n_steps = int(years*1/delta_t)
    
    #creating lists to store data
    x_list = [0.47] #list for x-values
    y_list = [0.0] #list for y-values
    v_x_list = [0.0] #list for x-velocities
    v_y_list = [8.17] #list for y- velocities
    ang_list = [] #list for angular momentum
    
    # setting constants
    M_sun = 1 # mass of sun in solar masses
    G = 39.5 # gravitational constant in solar masses, AU and yr as units
    
    #setting inital conditions for Mercury 
    x_initial_Merc = 0.47 #initial x-position of Mercury in AU 
    y_initial_Merc = 0.0 #initial y - position of Mercury in AU
    v_x_initial    = 0.0 #initial x-velocity in AU/yr
    v_y_initial    = 8.17 #initial y-veloctiy in AU/yr
    
    
    #loop to carry out numerical integration
    for i in range(1, n_steps+1):
     
     #defining 'r' 
     r_Merc = ((x_initial_Merc)**2 + (y_initial_Merc)**2)**0.5
      
     #new x velocity and position, uding Euler-Cromer integration 
     v_x_next = v_x_initial - ((G*M_sun*x_initial_Merc)/(r_Merc**3))*delta_t
     x_next  = x_initial_Merc + v_x_next*delta_t
    
     #new y velocity and position, uding Euler-Cromer integration 
     v_y_next = v_y_initial - ((G*M_sun*y_initial_Merc)/(r_Merc**3))*delta_t
     y_next   = y_initial_Merc + v_y_next*delta_t
     
     #finding angular momentum per unit mass, at each step, for Mercury 
     
     ang_mom_Merc = np.linalg.norm(np.cross((x_initial_Merc,y_initial_Merc),(v_x_initial,v_y_initial)))
     
     #adding data to previously created lists
     x_list.append(x_next) #list of x-values at each step
     y_list.append(y_next) #list of y-values at each step
     v_x_list.append(v_x_next) #list of x-velocities at each step
     v_y_list.append(v_y_next) #list of y-velocities at each step 
     ang_list.append(ang_mom_Merc)
     
     #re-setting values for next step of integration
     x_initial_Merc = x_next
     y_initial_Merc = y_next
     v_x_initial    = v_x_next
     v_y_initial    = v_y_next
     
    #checking if angular momentum is conserved
    if max(ang_list)-min(ang_list) != 0:
        print('angular momentum is not conserved but varies by:',max(ang_list)-min(ang_list),'times the mass of Mercury ' )
    else:
        print('angular momentum is conserved')
    #x vs y position plot for Mercury
    plt.figure()
    plt.errorbar(np.array([x_list]),np.array([y_list]),label='Mercury orbit')
    plt.plot(0,0, color = 'yellow', marker = 'o', label = "Sun", markersize = '5')
    plt.xlabel('x-position of Mercury in AU')
    plt.ylabel('y-position of Mercury in AU')
    plt.axis('equal')
    plt.legend(loc='upper right')
    plt.title('plot of x vs y position of Mercury over its orbit')
    plt.show()
    
    #x-component of velocity vs time graph for Mercury
    plt.figure()
    plt.errorbar(np.arange(0.0,years+delta_t,delta_t),v_x_list)
    plt.xlabel('Time in Earth year')
    plt.ylabel("x-component of Mercury's velocity in AU/yr" )
    plt.title('x-component of velocity vs time for Mercury')
    
    #y-component of velocity vs time graph for Mercury
    plt.figure()
    plt.errorbar(np.arange(0.0,years+delta_t,delta_t),v_y_list)
    plt.xlabel('Time in Earth year')
    plt.ylabel("y-component of Mercury's velocity in AU/yr" )
    plt.title('y-component of velocity vs time for Mercury')
    plt.show()

Mercury_orbit()

#1d
import matplotlib.pyplot as plt
import numpy as np
    
def Mercury_orbit():
    
    #setting time-step and number of steps for integration
    delta_t = 0.0001
    years = 1.0 #number of Earth years over which values are calculated 
    n_steps = int(years*(1/delta_t))
    
    #creating lists to store data
    x_list = [0.47] #list for x-values
    y_list = [0.0] #list for y-values
    v_x_list = [0.0] #list for x-velocities
    v_y_list = [8.17] #list for y- velocities
    ang_list = [] #list for angular momentum
    
    # setting constants
    M_sun = 1 # mass of sun in solar masses
    G = 39.5 # gravitational constant in solar masses, AU and yr as units
    alpha = 0.01 #relativisitic constant in AU^2
    #setting inital conditions for Mercury 
    x_initial_Merc = 0.47 #initial x-position of Mercury in AU 
    y_initial_Merc = 0.0 #initial y - position of Mercury in AU
    v_x_initial    = 0.0 #initial x-velocity in AU/yr
    v_y_initial    = 8.17 #initial y-veloctiy in AU/yr
    
    
    #loop to carry out numerical integration
    for i in range(1, n_steps+1):
     
     #defining 'r' for Mercury 
     r_Merc = ((x_initial_Merc)**2 + (y_initial_Merc)**2)**0.5
     
     #new x velocity and position, uding Euler-Cromer integration 
     v_x_next = v_x_initial - (((G*M_sun*x_initial_Merc)/(r_Merc**3))*(1 +(alpha/r_Merc**2))*delta_t)
     x_next  = x_initial_Merc + v_x_next*delta_t
    
     #new y velocity and position, uding Euler-Cromer integration 
     v_y_next = v_y_initial - (((G*M_sun*y_initial_Merc)/(r_Merc**3))*(1 +(alpha/r_Merc**2))*delta_t)
     y_next   = y_initial_Merc + v_y_next*delta_t
     
     #finding angular momentum per unit mass, at each step, for Mercury 
     ang_mom_Merc = np.linalg.norm(np.cross((x_initial_Merc,y_initial_Merc),(v_x_initial,v_y_initial)))
     
     #adding data to previously created lists
     x_list.append(x_next) #list of x-values at each step
     y_list.append(y_next) #list of y-values at each step
     v_x_list.append(v_x_next) #list of x-velocities at each step
     v_y_list.append(v_y_next) #list of y-velocities at each step 
     ang_list.append(ang_mom_Merc)
     
     #re-setting values for next step of integration
     x_initial_Merc = x_next
     y_initial_Merc = y_next
     v_x_initial    = v_x_next
     v_y_initial    = v_y_next
    
    #checking if angular momentum is conserved
    if max(ang_list)-min(ang_list) != 0:
        print('angular momentum is not conserved but varies by:',max(ang_list)-min(ang_list),'times the mass of Mercury ' )
    else:
        print('angular momentum is conserved')
    #x vs y position plot for Mercury
    plt.figure()
    plt.errorbar(np.array([x_list]),np.array([y_list]),label='Mercury orbit')
    plt.plot(0,0, color = 'yellow', marker = 'o', label = "Sun", markersize = '5')
    plt.xlabel('x-position of Mercury in AU')
    plt.ylabel('y-position of Mercury in AU')
    plt.axis('equal')
    plt.legend(loc='upper right')
    plt.title('plot of x vs y position of Mercury over its orbit\n(accounting for relativistic effects)')
    plt.show()
    
    #x-component of velocity vs time graph for Mercury
    plt.figure()
    plt.errorbar(np.arange(0.0,years+delta_t,delta_t),v_x_list)
    plt.xlabel('Time (in Earth year)')
    plt.ylabel("x-component of Mercury's velocity (in AU/yr)" )
    plt.title('x-component of velocity vs time for Mercury\n(accounting for relativistic effects)')
    
    #y-component of velocity vs time graph for Mercury
    plt.figure()
    plt.errorbar(np.arange(0.0,years+delta_t,delta_t),v_y_list)
    plt.xlabel('Time (in Earth year)')
    plt.ylabel("y-component of Mercury's velocity (in AU/yr)" )
    plt.title('y-component of velocity vs time for Mercury\n(accounting for relativistic effects)')
    plt.show()
    
Mercury_orbit()
