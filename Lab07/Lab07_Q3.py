from math import inf, sin,pi
import numpy as np#import empty,array,arange,isnan
import matplotlib.pyplot as plt

t_init = 0.0
t_final = 20.01
H = 20 # Size of inital step
delta = 1e-10     # Required position accuracy per unit time
n_max = 8
a = 1
b = 3
x = 0
y = 0
# recursion step
# else recall r t h, break it in half pg 3
def f(con):
    x = con[0]
    y = con[1]
    #print('f',x,y,1-(b+1)*x + a*x**2*y,(b)*x - a*x**2*y)
    vx = 1-(b+1)*x + a*x**2*y
    vy = b*x - a*x**2*y
    return np.array([vx,vy],float)

def step(r, t, H, n):
    r1 = r + H*f(r)

    # Do one modified midpoint step to get things started
    n = 1
    r1 = r + 0.5*H*f(r)
    r2 = r + H*f(r1)
    #print('r',r1,r2)

    # The array R1 stores the first row of the
    # extrapolation table, which contains only the single
    # modified midpoint estimate of the solution at the
    # end of the interval
    R1 = np.empty([1,2],float)
    R1[0] = 0.5*(r1 + r2 + 0.5*H*f(r2))

def full_calc(H):
    tpoints = np.arange(t_init,t_final,H)
    concentration = []
    r = np.array([x,y],float)

    # Do the "big steps" of size H

    for t in tpoints:

        concentration.append(r)
        #print(t,r)

        # Do one modified midpoint step to get things started
        n = 1
        r1 = r + 0.5*H*f(r)
        r2 = r + H*f(r1)
        #print('r',r1,r2)

        # The array R1 stores the first row of the
        # extrapolation table, which contains only the single
        # modified midpoint estimate of the solution at the
        # end of the interval
        R1 = np.empty([1,2],float)
        R1[0] = 0.5*(r1 + r2 + 0.5*H*f(r2))

        # Now increase n until the required accuracy is reached
        error = 2*H*delta
        while error>H*delta:
            
            n += 1
            if n > n_max:
                return full_calc(H/2)
            h = H/n
            #print("error",error,"h",h)

            # Modified midpoint method
            r1 = r + 0.5*h*f(r)
            r2 = r + h*f(r1)
            #print('in r',n,r1,r2)
            for i in range(n-1):
                r1 += h*f(r2)
                r2 += h*f(r1)
                #print('loop r',i,r1,r2)

            # Calculate extrapolation estimates.  Arrays R1 and R2
            # hold the two most recent lines of the table
            R2 = R1
            R1 = np.empty([n,2],float)

            R1[0] = 0.5*(r1 + r2 + 0.5*h*f(r2))
            # if isnan(R1[0][0]) or isnan(R1[0][1]):
            #     R1[0] = array([inf,inf])
            #print("EXCEPT",R1[0], r1, r2, h, f(r2))
            for m in range(1,n):
                epsilon = (R1[m-1]-R2[m-1])/((n/(n-1))**(2*m)-1)
                R1[m] = R1[m-1] + epsilon
            #print("eps",epsilon)
            if np.isnan(epsilon[0]) or np.isnan(epsilon[1]):
                epsilon = np.array([inf,inf])
            error = abs(epsilon[0])

        # Set r equal to the most accurate estimate we have,
        # before moving on to the next big step
        r = R1[n-1]
    #print(concentration)
    return tpoints, np.array(concentration)

tpoints, concentration = full_calc(H)
#print(tpoints,concentration)

# Plot the results
plt.figure()
plt.plot(tpoints,concentration[:,0],".-", label="x")
plt.plot(tpoints,concentration[:,1],".-", label="y")
plt.xlabel("Time (s)")
plt.ylabel("Chemical concentrations")
plt.title('Belousov-Zhabotinsky Reaction over Time')
plt.legend()
plt.grid()
plt.tight_layout()
plt.savefig('Q3.png', dpi=150)

plt.figure()
dtpoints = np.array(tpoints[1:]) - np.array(tpoints[:-1])
plt.plot(tpoints[:-1], dtpoints, label = "Timestep") # drop the last point in tpoints
plt.xlabel("Time t (s)")
plt.ylabel("Timestep h (s)")
plt.title('Size of Timestep over Time')
plt.legend()
plt.grid()
plt.tight_layout()

plt.show()
