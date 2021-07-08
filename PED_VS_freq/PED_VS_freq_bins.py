# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
from matplotlib import pyplot as plt

def calcbins(freqdata, binsize):
    mini = np.min(freqdata)
    maxi = np.max(freqdata)
    freqs = np.arange(np.round(mini,decimals=2),maxi+binsize,binsize)
    print freqs
    return freqs

filename = '2019-02-15.txt'
resonancefreq = 4252.72
binsize=0.004

freqdata, data = np.loadtxt(filename, usecols=(3,6),unpack=True)


actualdata = data>1 #for filtering out data points where cable was disconnected
data = data[actualdata]
freqdata = freqdata[actualdata]

#freqdata = np.around(freqdata, decimals=3) #for binning / filtering laser freq noise
#freqs = np.unique(freqdata)
freqs = calcbins(freqdata,binsize)
avgs = np.zeros(len(freqs))
stdevs = np.zeros(len(freqs))

for n in range(len(freqs)):
    cond1 = freqdata >= freqs[n]
    cond2 = freqdata < freqs[n]+binsize
    cond = cond1 * cond2
    avgs[n] = np.average(data[cond])
    stdevs[n] = np.std(data[cond])
    

plt.errorbar(freqs-resonancefreq, avgs, yerr=stdevs, fmt='o')
plt.title('CO2 beam energy VS excitation laser wavelength')
plt.xlabel('Laser wavelength deviation from resonance (nm)')
plt.ylabel('Lock-in voltage')
plt.show()
plt.close()
    