"""Created on 2021-09-14

@author: Gabriel Bailey

This code will manually sort N random samples into M linearly spaced bins, then do the same thing with numpy.histogram, and then displays a histogram of the results.
"""

import numpy as np
from time import time
import matplotlib.pyplot as plt

def histogram(min, max, M, samples):
    """Returns an array representing a histogram of val_array from min to max with M buckets.

    Args:
        min (float): The minimum value of the range
        max (float): The minimum value of the range
        M (integer): how many buckets to sort the values into
        val_array (array(float)): The array of data

    Returns:
        array(integer): The counts for each bucket
        array(integer): The edges of each bucket
    """
    #Initialize hist and bin_edges
    hist = np.zeros(M, dtype=int)
    bin_edges = np.arange(min,max,(max-min)/M, dtype=float)

    #Iterate accross all values in samples to count up hist
    for num in samples:
        i = int((num-min)/((max-min)/M)) #Use int() to truncate the value to the correct index
        if i >= 0 and i < M: #Ignore values outside the given bounds to prevent runtime errors
            hist[i] += 1
    return hist, bin_edges

#Set the parameters for the experiment
M = 1000
range_min = -5
range_max = 5
Ns = [10, 100, 1000, 10000, 100000, 1000000]
#The lists that store the times
times = []
nptimes = []

for N in Ns:
    #Create the sample data
    samples = np.random.randn(N)

    #Time how long it takes to run my histogram function
    start = time()
    hist, bin_edges = histogram(range_min, range_max, M, samples)
    end = time()

    #Time how long it takes to run numpy's histogram function
    npstart = time()
    nphist, npbin_edges = np.histogram(samples, bins=M, range=[range_min,range_max])
    npend = time()

    times.append(end-start)
    nptimes.append(npend-npstart)


    #Print the difference in time taken and whether the results are the same
    print("Time taken for", N, "samples (my hist vs np hist):", end-start, npend-npstart)
    print("histogram results are equal: ", np.array_equal(hist,nphist))

    #Plot both histogram results to compare (can comment this out if you want to test all 6 values of N immediately)
    """bar1 = plt.bar(bin_edges, hist, width=(range_max-range_min)/M, color = 'b',alpha=0.5)
    bar1.set_label("My histogram")
    bar2 = plt.bar(npbin_edges[:-1], nphist, width=(range_max-range_min)/M, color = 'r',alpha=0.5)
    bar2.set_label("numpy's histogram")
    plt.legend( title="Overlap is purple")
    plt.xlabel('Values')
    plt.ylabel('Frequency')
    plt.title('Histogram Comparision (N = '+str(N)+')')
    plt.yscale("log")
    plt.show()"""

# Plot a time comparision of our histogram function to numpy's histogram function with a logarithmic y axis.
plt.plot(Ns, times, color = 'b', label="My histogram")
plt.plot(Ns, nptimes, color = 'r', label = "numpy's histogram")
plt.legend()
plt.xlabel('N')
plt.ylabel('Time taken (s)')
plt.title('Histogram Function Time Comparision')
plt.yscale("log")
#plt.xscale("log")
plt.show()