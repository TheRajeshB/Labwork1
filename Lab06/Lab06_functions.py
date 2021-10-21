import numpy as np

#For question 2

# Finds the frequency of a list by assuming the start is a max and finding the next max
def find_freq_lazy(t,x):
    start = x[0]
    j = 100#len(x)//20
    while abs(x[j]) < abs(start*9999/10000):
        j += 1
    return 1/t[j]

# Finds the primary frequency of a list using fft
def find_freq_fft(maxtime,x):
    # Expand time domain by 100 times to increase frequency resolution:
    factor = 100
    c = np.fft.rfft(x,n=len(x)*factor)
    f = np.fft.rfftfreq(len(x)*factor,maxtime/len(x))
    maxi = np.argmax(abs(c))
    return f[maxi]
    # Can plot for sanity check:
    # fig = plt.figure(figsize=[10,5])
    # ax_freq = fig.add_subplot(1,1,1)
    # ax_freq.plot(f,abs(c))
    # ax_freq.set_xlabel('Frequency (Hz)')
    # ax_freq.set_ylabel('Magnitude ')

import colorsys
# From https://stackoverflow.com/a/49601444
# Just adjusts a color's lightness
def adjust_lightness(color, amount=0.5):
    import matplotlib.colors as mc
    import colorsys
    try:
        c = mc.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    return colorsys.hls_to_rgb(c[0], max(0, min(1, amount * c[1])), c[2])