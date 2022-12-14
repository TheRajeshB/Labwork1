'''
Authors: Gabriel Bailey and Farz Halim
This code simulates the position of Jupiter and Earth over time using Euler-Cromer integration and a Newtonian model. It then displays a trace of their positions, and a plot of the distance between them over time. Note that the Sun's position is not effected by Earth or Jupiter and Jupiter is not effected by Earth.

There are three cases this is done for:
a) 10 years simulation
b) 3 years simulation where Jupiter has the same mass as the Sun
c) 20 years simulation where Earth has the properties of an asteroid

See these cases being used at the bottom of the file.
'''
# Code uses the following libraries:
import matplotlib.pyplot as plt
import numpy as np

'''This functino simulates the system and plots the results for the given parameters. Note that duration is the only argument that is needed as all others have default values.
    '''
def Q2( duration,           # duration of simulation in years
        delta_t = 0.0001,   # timestep in seconds
        earth_x = 1.0,      # initial x-position of Earth in AU 
        earth_y = 0.0,      # initial y - position of Earth in AU
        earth_v_x = 0.0,    # initial x-velocity in AU/yr
        earth_v_y = 6.18,   # initial y-veloctiy in AU/yr
        jup_scale_factor = 1 # how much larger the mass of jupiter should be
        ):

    # setting constants
    M_sun = 1 # mass of sun in solar masses
    G = 39.5 # gravitational constant in solar masses, AU and yr as units
    M_jup = 0.001 * jup_scale_factor # mass of Jupiter in solar masses

    #setting time-step and number of steps for integration
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

    jup_r = ((jup_x)**2 + (jup_y)**2)**0.5 #distance from the Sun to Jupiter
    earth_r = ((earth_x)**2 + (earth_y)**2)**0.5 #distance from the Sun to Earth
    earth_jup_r = ((earth_x-jup_x)**2 + (earth_y-jup_y)**2)**0.5 #distance from Jupiter to Earth

    #loop to carry out numerical integration
    for i in range(0, n_steps):

        #adding data to storage lists
        times[i] = (i*delta_t)
        jup_x_list[i] = (jup_x) #list of x-values at each step
        jup_y_list[i] = (jup_y) #list of y-values at each step
        jup_v_x_list[i] = (jup_v_x) #list of x-velocities at each step
        jup_v_y_list[i] = (jup_v_y) #list of y-velocities at each step 
        jup_av_list[i] = (jup_av) #list of angular-velocities at each step

        earth_jup_r_list[i] = earth_jup_r
        earth_x_list[i] = (earth_x) #list of x-values at each step
        earth_y_list[i] = (earth_y) #list of y-values at each step
        earth_v_x_list[i] = (earth_v_x) #list of x-velocities at each step
        earth_v_y_list[i] = (earth_v_y) #list of y-velocities at each step 
        earth_av_list[i] = (earth_av) #list of angular-velocities at each step
        
        #Calculations Jupiter velocity 
        jup_v_x = jup_v_x - ((G*M_sun*jup_x)/(jup_r**3))*delta_t 
        jup_v_y = jup_v_y - ((G*M_sun*jup_y)/(jup_r**3))*delta_t
        
        #Calculations Earth's velocity based on gravitational attraction to the Sun and Jupiter
        earth_v_x = earth_v_x - ((G*M_sun*earth_x)/(earth_r**3) + (G*M_jup*(earth_x-jup_x))/(earth_jup_r**3))*delta_t 
        earth_v_y = earth_v_y - ((G*M_sun*earth_y)/(earth_r**3) + (G*M_jup*(earth_y-jup_y))/(earth_jup_r**3))*delta_t

        #finding angular momentum per unit mass, at each step, for Jupiter and Earth
        jup_av = np.linalg.norm(np.cross((jup_x,jup_y),(jup_v_x,jup_v_y))) 
        earth_av = np.linalg.norm(np.cross((earth_x,earth_y),(earth_v_x,earth_v_y)))

        #Calculations Jupiter's position, using Euler-Cromer integration
        jup_x   = jup_x + jup_v_x*delta_t
        jup_y   = jup_y + jup_v_y*delta_t

        #Calculations Earth's position, using Euler-Cromer integration
        earth_x   = earth_x + earth_v_x*delta_t
        earth_y   = earth_y + earth_v_y*delta_t
     

        #defining radiuses between the sun and the planets 
        jup_r = ((jup_x)**2 + (jup_y)**2)**0.5 #distance from the Sun to Jupiter
        earth_r = ((earth_x)**2 + (earth_y)**2)**0.5 #distance from the Sun to Earth
        earth_jup_r = ((earth_x-jup_x)**2 + (earth_y-jup_y)**2)**0.5 #distance from Jupiter to Earth
     
    #checking if angular momentum is conserved
    if max(earth_av_list)-min(earth_av_list) != 0:
        print('angular momentum is not conserved but varies by:',max(earth_av_list)-min(earth_av_list),'times the mass of Earth' )
    else:
        print('angular momentum is conserved')

    #x vs y position plot for Jupiter and Earth
    plt.figure()
    plt.plot(0,0, color = 'yellow', marker = 'o', label = "Sun", markersize = '5')
    plt.errorbar(jup_x_list, jup_y_list, color = 'orange', label = "Jupiter")
    plt.errorbar(earth_x_list, earth_y_list, color = 'b', label = "Earth")
    plt.legend()
    plt.xlabel('x-position in AU')
    plt.ylabel('y-position in AU')
    plt.axis('equal')
    plt.title('Plot of Position of Jupiter and Earth')  
    

    # Plot of distance from Earth to Jupiter over time
    plt.figure()
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
