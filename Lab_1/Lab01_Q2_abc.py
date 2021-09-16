'''
Authors: Gabriel Bailey and Farz Halim (?)
This code does:
'''

"""Created on Mon Sep 13 14:06:42 2021

@author: farzh
"""

import matplotlib.pyplot as plt
import numpy as np

def Q2( duration,
        delta_t = 0.0001, #timestep in seconds
        earth_x = 1.0, #initial x-position of Earth in AU 
        earth_y = 0.0, #initial y - position of Earth in AU
        earth_v_x = 0.0, #initial x-velocity in AU/yr
        earth_v_y = 6.18, #initial y-veloctiy in AU/yr
        jup_scale_factor = 1 #how much larger the mass of jupiter should be
        ):

    # setting constants
    M_sun = 1 # mass of sun in solar masses
    M_jup = 0.001 * jup_scale_factor # mass of Jupiter in solar masses
    G = 39.5 # gravitational constant in solar masses, AU and yr as units

    #setting time-step and number of steps for integration
    #duration = 10.0 # duration of simulation in Earth years
    #delta_t = 0.0001
    n_steps = int(duration/delta_t)

    #creating lists to store data
    times = np.empty(n_steps)
    earth_jup_r_list = np.empty(n_steps)
    jup_x_list = np.empty(n_steps) #list for x-values
    jup_y_list = np.empty(n_steps) #list for y-values
    jup_v_x_list = np.empty(n_steps) #list for x-velocities
    jup_v_y_list = np.empty(n_steps) #list for y- velocities
    jup_av_list = np.empty(n_steps) #list for angular momentum

    earth_x_list = np.empty(n_steps) #list for x-values
    earth_y_list = np.empty(n_steps) #list for y-values
    earth_v_x_list = np.empty(n_steps) #list for x-velocities
    earth_v_y_list = np.empty(n_steps) #list for y- velocities
    earth_av_list = np.empty(n_steps) #list for angular momentum
    
    #setting inital conditions for Jupiter 
    jup_x = 5.2 #initial x-position of Jupiter in AU 
    jup_y = 0.0 #initial y - position of Jupiter in AU
    jup_v_x = 0.0 #initial x-velocity in AU/yr
    jup_v_y = 2.63 #initial y-veloctiy in AU/yr
    jup_av = np.linalg.norm(np.cross((jup_x,jup_y),(jup_v_x,jup_v_y)))

    #setting inital conditions for Earth
    #earth_x = 1.0 #initial x-position of Earth in AU 
    #earth_y = 0.0 #initial y - position of Earth in AU
    #earth_v_x = 0.0 #initial x-velocity in AU/yr
    #earth_v_y = 6.18 #initial y-veloctiy in AU/yr
    earth_av = np.linalg.norm(np.cross((earth_x,earth_y),(earth_v_x,earth_v_y))) #initial angular momentum in AU^2/yr

    #loop to carry out numerical integration
    for i in range(0, n_steps):

        #adding data to previously created lists
        times[i] = (i*delta_t)
        jup_x_list[i] = (jup_x) #list of x-values at each step
        jup_y_list[i] = (jup_y) #list of y-values at each step
        jup_v_x_list[i] = (jup_v_x) #list of x-velocities at each step
        jup_v_y_list[i] = (jup_v_y) #list of y-velocities at each step 
        jup_av_list[i] = (jup_av) #list of angular-velocities at each step

        #defining 'r' 
        jup_r = ((jup_x)**2 + (jup_y)**2)**0.5 #distance from the Sun to Jupiter
        earth_r = ((earth_x)**2 + (earth_y)**2)**0.5 #distance from the Sun to Earth
        earth_jup_r = ((earth_x-jup_x)**2 + (earth_y-jup_y)**2)**0.5 #distance from Jupiter to Earth
        
        #new x velocity and position, uding Euler-Cromer integration 
        jup_v_x = jup_v_x - ((G*M_sun*jup_x)/(jup_r**3))*delta_t
        jup_x   = jup_x + jup_v_x*delta_t

        #new y velocity and position, uding Euler-Cromer integration 
        jup_v_y = jup_v_y - ((G*M_sun*jup_y)/(jup_r**3))*delta_t
        jup_y   = jup_y + jup_v_y*delta_t

        #finding angular momentum per unit mass, at each step, for Jupiter 
        jup_av = np.linalg.norm(np.cross((jup_x,jup_y),(jup_v_x,jup_v_y)))

        #Calculations for Earth
        earth_jup_r_list[i] = earth_jup_r
        earth_x_list[i] = (earth_x) #list of x-values at each step
        earth_y_list[i] = (earth_y) #list of y-values at each step
        earth_v_x_list[i] = (earth_v_x) #list of x-velocities at each step
        earth_v_y_list[i] = (earth_v_y) #list of y-velocities at each step 
        earth_av_list[i] = (earth_av) #list of angular-velocities at each step

        #new x velocity and position, uding Euler-Cromer integration 
        earth_v_x = earth_v_x - ((G*M_sun*earth_x)/(earth_r**3) + (G*M_jup*(earth_x-jup_x))/(earth_jup_r**3))*delta_t
        earth_x   = earth_x + earth_v_x*delta_t

        #new y velocity and position, uding Euler-Cromer integration 
        earth_v_y = earth_v_y - ((G*M_sun*earth_y)/(earth_r**3) + (G*M_jup*(earth_y-jup_y))/(earth_jup_r**3))*delta_t
        earth_y   = earth_y + earth_v_y*delta_t

        #finding angular momentum per unit mass, at each step, for Earth 
        earth_av = np.linalg.norm(np.cross((earth_x,earth_y),(earth_v_x,earth_v_y)))
     
     
    #checking if angular momentum is conserved
    if max(earth_av_list)-min(earth_av_list) != 0:
        print('angular momentum is not conserved but varies by:',max(earth_av_list)-min(earth_av_list),'times the mass of Earth' )
    else:
        print('angular momentum is conserved')

    #x vs y position plot for Jupiter and Earth
    plt.plot(0,0, color = 'yellow', marker = 'o', label = "Sun", markersize = '5')
    plt.errorbar(jup_x_list, jup_y_list, color = 'orange', label = "Jupiter")
    plt.errorbar(earth_x_list, earth_y_list, color = 'b', label = "Earth")
    plt.legend()
    plt.xlabel('x-position in AU')
    plt.ylabel('y-position in AU')
    plt.title('Plot of Position of Jupiter and Earth')
    plt.show()

    # Plot of distance from Earth to Jupiter over time
    plt.errorbar(times,earth_jup_r_list, color = 'b')
    plt.xlabel('Time in Earth years')
    plt.ylabel('Distance in AU')
    plt.title('Distance Between Jupiter and Earth Over Time')
    plt.show()

#2a
Q2(duration = 10.0, # duration of simulation in years
    )

#2b
Q2(duration = 3.0, # duration of simulation in Earth years
    jup_scale_factor = 1000
    )

#2c
Q2(duration = 20.0, # duration of simulation in Earth years
    earth_x = 3.3, #initial x-position of Earth in AU 
    earth_y = 0.0, #initial y - position of Earth in AU
    earth_v_x = 0.0, #initial x-velocity in AU/yr
    earth_v_y = 3.46, #initial y-veloctiy in AU/yr
    )


#FARZ'S VERSION 

import numpy as np 
import matplotlib.pyplot as plt 

def three_body_problem():
    
    #setting up time step for integration 
    delta_t = 0.0001
    years   = 3.0 #number of years over which values are calculated
    n_steps = int(years*1/delta_t) #number of steps for numerical integration
    
    #creating lists to store output at each step of EC integration 
    E_position_x = [1.0] #list to store x-position of Earth
    E_position_y = [0.0] #list to store y-position of Earth
    J_position_x = [5.2] #list to store x-position of Jupiter
    J_position_y = [0.0] #list to store y-position of Jupiter 
    
    
    
    #setting constants 
    M_sun = 1 # mass of sun in solar masses
    G = 39.5 # gravitational constant in solar masses, AU and yr as units
    M_Jupiter = 1 #mass of Jupiter in solar masses
    
    #setting initial conditions for Earth 
    x_E_initial     =  1.0 #initial Earth x-position in AU 
    y_E_initial     =  0.0 #initial Earth y-psoition in AU 
    v_E_initial_x   =  0.0 #initial x-component of Earth velocity in AU/yr
    v_E_initial_y   =  6.18 #initial y-component of Earth velocity in AU/yr
    
    #setting initial conditions for Jupiter
    x_J_initial   = 5.2  #initial Jupiter x-position in AU 
    y_J_initial   = 0.0  #initial Jupiter y-position in AU
    v_J_initial_x = 0.0  #initial x-component of Jupiter velocity in AU/yr
    v_J_initial_y = 2.63 #initial y-component of Jupiter velocity in AU/yr
    
    #loop to carry out numerical integration
    for i in range(0,n_steps):
        
        #defining distance of Jupiter from Sun(r_J)
        r_J = ((x_J_initial)**2 + (y_J_initial)**2)**0.5
        
        #defining distance of Earth from Sun (r_E)
        r_E = ((x_E_initial)**2 + (y_E_initial)**2)**0.5 
        
        #defining distance between Earth and Jupiter
        r_E_J = ((x_E_initial-x_J_initial)**2 + (y_E_initial-y_J_initial)**2)**0.5
        
        #determining new velocity components for Earth using EC method
        v_E_next_x = v_E_initial_x - (G*M_sun*x_E_initial/r_E**3)*delta_t + (G*M_Jupiter*(x_J_initial-x_E_initial)/r_E_J**3)*delta_t
        v_E_next_y = v_E_initial_y - (G*M_sun*y_E_initial/r_E**3)*delta_t + (G*M_Jupiter*(x_J_initial-x_E_initial)/r_E_J**3)*delta_t
        
        #determining new velocity components for Jupiter, ignoring influence of Earth 
        v_J_next_x = v_J_initial_x - (G*M_sun*x_J_initial)/r_J**3*delta_t
        v_J_next_y = v_J_initial_y - (G*M_sun*y_J_initial)/r_J**3*delta_t
        
        #determining new components for Earth position
        x_E_next = x_E_initial + v_E_next_x*delta_t
        y_E_next = y_E_initial + v_E_next_y*delta_t
        
        #determining new components for Jupiter position
        x_J_next = x_J_initial + v_J_next_x*delta_t
        y_J_next = y_J_initial + v_J_next_y*delta_t
        
        #storing outputs in previously created lists
        E_position_x.append(x_E_next)
        E_position_y.append(y_E_next)
        J_position_x.append(x_J_next)
        J_position_y.append(y_J_next)
        
        #re-setting values for next step of integration
        x_E_initial     = x_E_next
        y_E_initial      = y_E_next
        v_E_initial_x   = v_E_next_x
        v_E_initial_y   = v_E_next_y
        
        x_J_initial     = x_J_next
        y_J_initial     = y_J_next
        v_J_initial_x   = v_J_next_x
        v_J_initial_y   = v_J_next_y
    
    plt.errorbar(E_position_x,E_position_y)
    plt.errorbar(J_position_x,J_position_y)
    plt.show()

print(three_body_problem())
