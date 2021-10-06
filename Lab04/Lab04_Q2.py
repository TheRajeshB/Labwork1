"""Created on 2021-10-06

@author: Gabriel Bailey

This code will calculate the eigenvalues for an electron in an energy well V(x) = a*x/L with infinitely high walls, and print out the first 10 energy levels, and plot the wavefunctions for the first 3 energy levels.
"""

import numpy as np
from scipy.constants import pi,electron_volt,physical_constants,electron_mass
from functions_Lab03 import simp_integrate
import matplotlib.pyplot as plt

# Define constants
hbar = physical_constants["Planck constant over 2 pi in eV s"][0]

L = 5 * 10**-10 #m
a = 10 #eV
#M = 9.1094 * 10**-31 #kg
M = electron_mass
q = 1.6022 * 10**-19 #C
print(M)

#Define the values of the H matrix
def Hmatrix(m, n):
    if m != n:
        if m+n % 2 == 0:
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

# Calculate the eigenvalues (Energies in eV)
eigen_vals10 = np.linalg.eigvals(H10)
eigen_vals100 = np.linalg.eigvals(H100)

# Print the results
print("Eigenvalues 10 :",np.sort(eigen_vals10))
print("Eigenvalues 100:", np.sort(eigen_vals100)[0:10])

eigvals, eigvecs = np.linalg.eig(H100)
eigorder = np.argsort(eigvals) #Fins the sorted order from smallest to largest

print("Eigenvalues:", eigvals[eigorder][0:3])
#print("Eigenvectors:", eigvecs[:,eigorder][0:3])

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
        psi += psi_n[n] * np.sin(pi*n*x/L)
    return psi

# Define domains (101 so that we include L for simp integration)
x = np.linspace(0,L,101,True)
firstx = 3

# Calculate wavefuntions
wave = np.empty((firstx,N+1))
for i in range(firstx):
    # Calculate values
    wave[i] = psi(eigvecs[:,eigorder][i],x)
    #Normalize
    wave[i] /= np.sqrt(simp_integrate((lambda y : abs(wave[i][int(N*y//L)])**2), N, 0, L))

V = a*x/L*np.average(wave[0])/8
# Plot the results
fig1 = plt.figure(figsize=[10,5])
ax = fig1.add_subplot(1,1,1)
ax.plot(x, wave[0], '-', label = "N = 1")
ax.plot(x, wave[1], '-', label = "N = 2")
ax.plot(x, wave[2], '-', label = "N = 3")
ax.plot(x, V, '--', label = "V (scaled)")
ax.legend()
ax.grid()
ax.set_xlabel('x (m)')
ax.set_ylabel('Wavefunction Value')
ax.set_title('Energy states for electron with V=ax/L')  

plt.show()

