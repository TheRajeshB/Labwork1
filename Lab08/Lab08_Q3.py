# Started from laplace.py

from numpy import arange, empty, pi,zeros,max, sin,linspace
from pylab import imshow,gray,show
import matplotlib.pyplot as plt

# Constants
epsilon = 1      # Target accuracy
delta_x = 0.02      # 
delta_t = 0.005     #
Lx = 2*pi
Tf = 3

Nx = int(Lx//delta_x)
Nt = int(Tf//delta_t)
# Create arrays to hold potential values
u = zeros([Nx, Nt],float)
# Set initial conditions
u[0,:] = 0
u[Nx-1,:] = 0
for x in range(1,Nx-1):
    u[x,0] = sin(x*delta_x)
    # Also forward Euler step:
    u[x,1] = u[x,0] - epsilon*delta_t/(4*delta_x)*(u[x+1,0]**2-u[x-1,0]**2)

beta = epsilon*delta_t/delta_x
for t in range(1,Nt-1):
    for x in range(1,Nx-1):
        u[x,t+1] = u[x,t-1] - beta/2*(u[x+1,t]**2-u[x-1,t]**2)
        '''try:
            u[x,t+1] = u[x,t-1] - beta/2*(u[x+1,t]**2-u[x-1,t]**2)
        except:
            print(t,x,u[x,t-1],u[x+1,t],u[x-1,t])'''

# at t = 0, 0.5, 1, 1.5
# Make a plot
#print(Nx,Nt)
def plot_result(u,time):
    t = int(time/Tf * Nt)
    #print(t)
    fig = plt.figure(figsize=[10,5])
    ax = fig.add_subplot(1,1,1)
    ax.set_xlim(0, Lx)
    ax.set_ylim(-1.6, 1.6)
    x = linspace(0, Lx, Nx)
    ax.plot(x,u[:,t])
        
    ax.grid(which='both', axis='y')
    ax.set_title('Burger\'s equation at Time {}s'.format(time)) 
    ax.set_xlabel('x Space')
    ax.set_ylabel('Value')
    return ax

plot_result(u,0.0)
plot_result(u,0.5)
plot_result(u,1.0)
plot_result(u,1.5)
#plt.show()
# Plot it as an animation:

from matplotlib.animation import FuncAnimation

fig, ax = plt.subplots()
ax.grid(which='both', axis='y')
ax.set_title('Burger\'s equation') 
ax.set_xlabel('x Space')
ax.set_ylabel('Value')
xdata, ydata = [], []
line, = plt.plot([], [])
time_text = ax.text(0.1, 0.1, 'time = ', transform=ax.transAxes)

def init():
    ax.set_xlim(0, 2*pi)
    ax.set_ylim(-2, 2)
    time_text.set_text('time = ')
    return line,

def update(i):
    i = int(i)
    x = linspace(0, Lx, Nx)
    #print(i)
    y = u[:,i]
    line.set_data(x, y)
    #print(i*delta_t)
    time_text.set_text('time = {:1.2f}'.format(i*delta_t))
    return line, time_text

ani = FuncAnimation(fig, update, frames=Nt, interval=1,
                    init_func=init, blit=True)
plt.show()
