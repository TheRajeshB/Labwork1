'''Created on 2021-10-20

@author: Gabriel Bailey

This code will simulate the motion of a building using the verlet method plotting the result, and then calculate the normal modes and frequencies.
'''
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import colorsys

def generateA(N):
    A = np.zeros((N,N))
    for i in range(N):
        A[i][i] = -2
        if i > 0:
            A[i][i-1] = 1
        if i < N-1:
            A[i][i+1] = 1
    return A

def simulate_building(N,maxtime,h,x0):
    steps = int(maxtime/h)
    time = np.linspace(0,maxtime,steps)
    x = np.zeros((2*steps,N))
    x[0][0] = x0
    v = np.zeros((2*steps,N))
    A = km*generateA(N)
    #first step
    v[1] = v[0]+1/2*h*A.dot(x[0])

    for t in range(2,2*steps,2):
        x[t] = x[t-2]+h*v[t-1]
        k = h*A.dot(x[t])
        v[t] = v[t-1]+1/2*k
        v[t+1] = v[t-1] + k
        # if t < 6:
        #     print(x[t])
        #     print(k)
        #     print(v[t])

    return time, x, v

# Define all the constants for this problem (all in m)
km = 400 #rad s^-2
dt = 1/1000 #s
maxtime = 1 #s
x0 = 10 # cm
N = 10
t10,x10,v10 = simulate_building(N,maxtime,dt,x0)

#print(x10[:20:2,:].T)

# From https://stackoverflow.com/a/49601444
# Just adjusts a colors lightness
def adjust_lightness(color, amount=0.5):
    import matplotlib.colors as mc
    import colorsys
    try:
        c = mc.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    return colorsys.hls_to_rgb(c[0], max(0, min(1, amount * c[1])), c[2])

# Plot the results
fig = plt.figure(figsize=[10,5])
ax = fig.add_subplot(1,1,1)
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
for i in range(N):
    ax.axhline(y=10*i, linestyle=':', color = adjust_lightness(colors[i],0.7))
    ax.plot(t10, x10[::2,i].T+10*i)
    
#ax.set_xlim( -0.05, 0.05)
ax.grid(which='both', axis='y')
ml = MultipleLocator(10)
ax.yaxis.set_minor_locator(ml)
#ax.tick_params(axis='x', which='minor', bottom=False)
ax.set_title('Vibrations of a Building with {} Floors\n (dotted line is 0 for each floor)'.format(N)) 
ax.set_xlabel('Time (s)')
ax.set_ylabel('x displacement for floor 0 (cm) ')
    
plt.show()

#Part b

