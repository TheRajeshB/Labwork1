# -*- coding: utf-8 -*-
'''Created on 2021-11-23

@author: Gabriel Bailey

This code will simulate photons leaving the photosphere of a star to generate the angle at which they leave and the limb-darkening effect. It will produce 2 plot: One for the final scattering angle and one for the intensity. The variable 'part' determines whether the results for part b or part c are generated.
'''
from logging import error
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

# The function provided:
def get_tau_step():
    """ Calculate how far a photon travels before it gets scattered.
    OUT: optical depth traveled """
    delta_tau = -np.log(np.random.random())
    return delta_tau

def emit_photon(tau_max):
    """ Emit a photon from the stellar core.
    IN: tau max is max optical depth
    OUT:
    tau: optical depth at which the photon is created
    mu: directional cosine of the photon emitted """
    tau = tau_max
    delta_tau = get_tau_step()
    mu = np.random.random()
    return tau - delta_tau*mu, mu

def scatter_photon(tau):
    """ Scatter a photon.
    IN: tau, optical depth of the atmosphere
    OUT:
    tau: new optical depth
    mu: directional cosine of the photon scattered """
    delta_tau = get_tau_step()
    mu = 2*np.random.random()-1 # sample mu uniformly from -1 to 1
    tau = tau - delta_tau*mu
    return tau, mu

# Culmative function that calculate the final mu of a photon leaving the star
def full_scatter(tau_max):
    tau, mu = emit_photon(tau_max)
    count = 0
    while tau >= 0:
        if tau > tau_max: # reset to 0 if it goes back to the core
            return full_scatter(tau_max) #redo it
        tau, mu = scatter_photon(tau)
        count += 1
    return mu, count

# options
compute_scatters = True
part = 'c' # Choose between part b and c
N = int(1e5)

tau_max = 10
if part == 'b':
    tau_max = 10
elif part == 'c':
    tau_max = 1e-4
else:
    error('Invalid part of question selected (must be \'b\' or \'c\').')

# Array to put data into
mus = np.empty(N)

# Save data rather than recalculate it if desired:
if compute_scatters:
    for i in range(N):
        mus[i] ,count = full_scatter(tau_max)

    np.savez('Q2_mus_'+part, mus=mus)
else:
    npzfile = np.load('Q2_mus_'+part+'.npz')
    mus = npzfile['mus']

# For testing if results are in the expected bounds
""" print(np.min(mus),np.max(mus))
for i in range(len(mus)):
    if mus[i]<0:
        print(mus[i]) """

# Calculate the properties of the mu histogram
hist, bin_edges = np.histogram(mus, bins = 20, range=(0,1))
Ns = hist.copy()
hist = hist/np.max(hist)
bin_loc = (bin_edges[:-1]+bin_edges[1:])/2
width = 0.8*(bin_edges[1]-bin_edges[0])

# Plot the resulting mu histogram
fig, axb1 = plt.subplots()
axb1.set_title('Final Scattering Angle') 
axb1.set_xlabel('Mu')
axb1.set_ylabel('N (scaled)') 
axb1.bar(bin_loc, hist, width, align = 'center')

# Calculate the intensity
Is = Ns/bin_loc
Is = Is/np.max(Is)

# Curvefit the intensity
def calc_I(mu,a,b):
    return (a+b*mu)
popt, pcov = curve_fit(calc_I, bin_loc, Is)

# Plot the results for intensity
fig, axb2 = plt.subplots()
axb2.set_title('Limb-Darkening Intensity')
axb2.set_xlabel('Mu')
axb2.set_ylabel('Intensity (Scaled)')
axb2.bar(bin_loc, Is, width, align = 'center', label='Data Histogram')
axb2.plot(bin_loc,calc_I(bin_loc,popt[0],popt[1]),label='Best Fit',color='red')
axb2.plot(bin_loc,calc_I(bin_loc,0.4,0.6),label='Analytic Fit', color='orange')
axb2.text(0.05,0.7,'Best fit:\nI=({:.2f}+{:.2f}*mu)/I_1'.format(popt[0],popt[1]))
axb2.legend()

plt.show()