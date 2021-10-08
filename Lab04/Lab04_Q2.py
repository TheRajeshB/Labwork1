"""Created on 2021-10-06

@author: Gabriel Bailey

This code will calculate the eigenvalues for an electron in an energy well V(x) = a*x/L with infinitely high walls, print out the first 10 energy levels, and plot the wavefunctions for the first 3 energy levels.
"""

import numpy as np
from scipy.constants import pi,electron_volt,physical_constants,electron_mass
from Lab04_functions import simp_integrate
import matplotlib.pyplot as plt

# Define constants
hbar = physical_constants["Planck constant over 2 pi in eV s"][0]

L = 5 * 10**-10 #m
a = 10 #eV

#M = 9.1094 * 10**-31 #kg, textbook definition
M = electron_mass # 9.1093837015e-31
q = 1.6022 * 10**-19 # Charge

#Define the values of the H matrix
def Hmatrix(m, n):
    if m != n:
        if (m+n) % 2 == 0:
            return 0
        else:
            return - (8*a*m*n)/(pi**2 * (m**2-n**2)**2)
    else:
        return 1/2*a + pi**2 * hbar**2 * m**2/(2*M*L**2) * electron_volt

# Calculate the H matrix
def calc_H(mmax,nmax):
    H = np.empty((mmax,nmax))
    for m in range(1, mmax+1):
        for n in range(1, nmax+1):
            H[m-1, n-1] = Hmatrix(m, n)
    return H

# Calculate the values of H for different Ns
N = 10
H10 = calc_H(N,N)
N = 100
H100 = calc_H(N,N)
'''
print("Part c and d")
# Calculate the eigenvalues (Energies in eV)
eigen_vals10 = np.linalg.eigvals(H10)
eigen_vals100 = np.linalg.eigvals(H100)

# Print the results
print("Eigenvalues for N = 10 :",np.sort(eigen_vals10))
print("Eigenvalues for N = 100:", np.sort(eigen_vals100)[0:10])

# Estimating proportional errors
#print(','.join(["{:.4f}".format(fs) for fs in np.log10(abs((np.sort(eigen_vals10)-np.sort(eigen_vals100)[:10])/np.sort(eigen_vals10)))]))
'''

print("Part e")
# Calculate the eigenvalues and eigenvectors
eigvals, eigvecs = np.linalg.eig(H100)
# np docs: "the column v[:,i] is the eigenvector corresponding to the eigenvalue w[i]"

# Order them in terms of energy (from the eigenvalues)
eigorder = np.argsort(eigvals) # An array of indices such that eigvals[eigorder] yields a sorted eigvals
eigvals = eigvals[eigorder]
eigvecs = eigvecs[:,eigorder] # Apply it to the column order of eigvecs

firstx = 3
print("Eigenvalues:", eigvals[:firstx])
#print("Eigenvectors:", eigvecs[:firstx])

def psi(psi_n,x):
    """Calculate the value of the wavefunction at x given some Psi_n

    Args:
        psi_n (numpy.Array): An array of the values of psi_n
        x (float): A position to calculate the wavefunction at

    Returns:
        float: The value of the wavefunction at x
    """
    psi = 0
    for n in range(len(psi_n)):
        psi += psi_n[n] * np.sin(pi*(n+1)*x/L)
    return psi

# Define domains (101 so that we include L for simp integration)
x = np.linspace(0,L,101,True)

# Calculate wavefuntions
wave = np.empty((firstx,N+1))
for i in range(firstx):
    # Calculate values
    wave[i] = psi(eigvecs[:,i],x)
    #Normalize
    A = simp_integrate((lambda y : abs(wave[i][int(N*y//L)])**2), N, 0, L)
    wave[i] /= np.sqrt(A)

# Trash method to get i from y easily
def ytoi(y, start, stop):
    for i in range(len(x)):
        if x[i] >= start and x[i] <= stop:
            if abs(x[i]-y) <= L/N/2:
                return i
    return -1

# Create a scaled V function so we can see the well
V = np.array([0,0,a*np.max(wave)/20,0])
V[[0,-1]] = np.max(wave)
# x-values to go with it
Vx = np.array([0,0+L/1000,L-L/1000,L])

# Plot the results
fig1 = plt.figure(figsize=[10,5])
ax = fig1.add_subplot(1,1,1)
for i in range(firstx):
    ax.plot(x, wave[i], '-', label = "n = "+str(i+1)) #+", E = "+"{:.2f}".format(eigvals[i])) # If we want energy on plot too
ax.plot(Vx, V, '--', label = "V (scaled)")
ax.legend(title="Energy Levels")
ax.grid()
ax.set_xlabel('x (m)')
ax.set_ylabel('Wavefunction Value')
ax.set_title('Energy States for Electron with V=ax/L')  

plt.show()

