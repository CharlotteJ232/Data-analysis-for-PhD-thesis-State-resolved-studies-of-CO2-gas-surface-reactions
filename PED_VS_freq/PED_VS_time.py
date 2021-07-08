# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
from matplotlib import pyplot as plt

filename = '2019-02-15'
ext = '.txt'
resonancefreq = 4252.72

plotfreq = 4252.724
freqvariation = 0.002

datalowV = 1
datahighV = 100


freqdata, data = np.loadtxt('Data/'+filename+ext, usecols=(3,6),unpack=True)

actual1 = data>datalowV #for filtering out data points where cable was disconnected
actual2 = data<datahighV
actualdata = actual1*actual2
data = data[actualdata]
freqdata = freqdata[actualdata]

cond1 = freqdata > plotfreq-freqvariation
cond2 = freqdata < plotfreq+freqvariation
cond = cond1*cond2
    

plt.plot(np.arange(len(data[cond])),data[cond])
plt.title('CO2 beam energy VS time for fixed wavelength: '+
          str(plotfreq)+'+-'+str(freqvariation)+' nm')
plt.xlabel('Time (not continuous and arbitrary units)')
plt.ylabel('Lock-in voltage')
plt.savefig('Figures/'+str(filename)+'_'+str(plotfreq)+'+-'+str(freqvariation)+'.png',dpi=1000)
plt.show()
plt.close()
    