# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
from matplotlib import pyplot as plt

filename = '2019-02-15.txt'
resonancefreq = 4252.72
decimal = 3
datalowV = 1
datahighV = 4

freqdata, data = np.loadtxt(filename, usecols=(3,6),unpack=True)

actual1 = data>datalowV #for filtering out data points where cable was disconnected
actual2 = data<datahighV
actualdata = actual1*actual2
data = data[actualdata]
freqdata = freqdata[actualdata]

freqdata = np.around(freqdata, decimals=decimal) #for binning / filtering laser freq noise

freqs = np.unique(freqdata)
avgs = np.zeros(len(freqs))
stdevs = np.zeros(len(freqs))

for n in range(len(freqs)):
    cond = freqdata == freqs[n]
    
    avgs[n] = np.average(data[cond])
    stdevs[n] = np.std(data[cond])
    

plt.errorbar(freqs-resonancefreq, avgs, yerr=stdevs, fmt='o')
plt.title('CO2 beam energy VS excitation laser wavelength')
plt.xlabel('Laser wavelength deviation from resonance (nm)')
plt.ylabel('Lock-in voltage')
plt.show()
plt.close()
    