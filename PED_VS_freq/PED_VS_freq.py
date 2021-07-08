# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
from __future__ import division
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import cm

filename = '2019-02-15'
ext = '.txt'
resonancefreq = 4252.72
assen = [-0.01, 0.01, None,None]
datalowV = 5
datahighV = 40

freqdata, data = np.loadtxt('Data/'+filename+ext, usecols=(3,6),unpack=True)

actual1 = data>datalowV #for filtering out data points where cable was disconnected
actual2 = data<datahighV
actualdata = actual1*actual2
data = data[actualdata]
freqdata = freqdata[actualdata]
   
col = cm.get_cmap('rainbow', len(data)) 
collist = [col(i) for i in range(col.N)]
plt.scatter(freqdata-resonancefreq, data, marker='.',color=collist)
plt.axis(assen)
plt.title('CO2 beam energy VS excitation laser wavelength')
plt.xlabel('Laser wavelength deviation from resonance (nm)')
plt.ylabel('Lock-in voltage')
plt.savefig('Figures/'+str(filename)+'_'+str(assen)+'_'+str(datalowV)+'-'+str(datahighV)+'.png',dpi=1000)
plt.show()
plt.close()
    