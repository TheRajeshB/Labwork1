'''
Authors: Gabriel Bailey and Farz Halim (?)
This code will histogram N random samples into M linearly spaced bins
'''
import numpy as np
from time import time
import matplotlib.pyplot as plt

def histogram(min, max, M, val_array):
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
    hist = np.zeros(M, dtype=int)
    bin_edges = np.arange(min,max,(max-min)/M, dtype=float)
    for num in val_array:
        i = int((num-min)/((max-min)/M))
        if i >= 0 and i < M:
            #print(i,(max-min)/M,num-min, (num-min)/((max-min)/M))
            hist[i] += 1
    return hist, bin_edges

M = 1000
min = -5
max = 5
for N in [10, 100, 1000, 10000, 100000, 1000000]:
    samples = np.random.randn(N)

    start = time()
    hist, bin_edges = histogram(min, max, M, samples)
    end = time()

    npstart = time()
    nphist, npbin_edges = np.histogram(samples, bins=1000, range=[min,max])
    npend = time()

    #print(end-start, npend-npstart)
    print("Time taken for", N, "samples (my hist vs np hist):", end-start, npend-npstart)
    print("histogram results are equal: ", np.array_equal(hist,nphist))
    bar1 = plt.bar(bin_edges, hist, width=(max-min)/M, color = 'b',alpha=0.5)
    bar1.set_label("My histogram")
    bar2 = plt.bar(npbin_edges[:-1], nphist, width=(max-min)/M, color = 'r',alpha=0.5)
    bar2.set_label("numpy's histogram")
    plt.legend()
    plt.xlabel('Values')
    plt.ylabel('Frequency')
    plt.title('Histogram Comparision (overlap is purple)')
    plt.yscale("log")
    plt.show()
'''
Results:
0.0
0.000965118408203125
0.0019991397857666016
0.025830745697021484
0.21095919609069824
2.0834832191467285
'''
'''Results for my hist vs np hist:
0.0 0.0
0.0 0.0
0.0 0.0038881301879882812
0.009621143341064453 0.0010001659393310547
0.08453035354614258 0.0010955333709716797
0.8508334159851074 0.01567864418029785
My method is slower, likely due to np's methods directly using C code while I use numpy's methods to access arrays. My method seems to have slightly worse time complexity: the proportional difference between them grows with n (from ~10x faster at 10k to ~50x faster at 1 million). This is strange since my method is O(N).'''

