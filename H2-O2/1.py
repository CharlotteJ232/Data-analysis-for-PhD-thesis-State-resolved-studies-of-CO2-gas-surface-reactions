# -*- coding: utf-8 -*-
"""
Created on Thu Feb 21 14:15:02 2019

@author: Charlotte
"""

from __future__ import division
import numpy as np
from matplotlib import pyplot as plt

file1 = 'h2_co2'
file2 = 'o2_co2'
file3 = 'ar_co2_co'
file5 = 'ar_o2_co2_co_cont'
file4 = 'ar_o2_co2_co_disc'
ext = '.asc'

time1, data1 = np.loadtxt('Data/'+file1+ext, usecols=(4,5), 
                                skiprows=14, unpack=True)                                                                               
time2, data2 = np.loadtxt('Data/'+file2+ext, usecols=(4,5), 
                                skiprows=14, unpack=True) 
time3, data3 = np.loadtxt('Data/'+file3+ext, usecols=(4,5), 
                                skiprows=15, unpack=True) 
time4, data4 = np.loadtxt('Data/'+file4+ext, usecols=(4,5), 
                                skiprows=16, unpack=True) 
time5, data5 = np.loadtxt('Data/'+file5+ext, usecols=(4,5), 
                                skiprows=16, unpack=True)                                 


#-----------------------------------------------------------------------------#                
                               
for time, data, lab in zip([time1, time2, time3, time4, time5], 
                             [data1, data2, data3, data4, data5], 
                             [file1, file2, file3, file4, file5]):
    data = data/np.max(data[0:200])
    plt.plot(time/60, data, label=lab)
plt.legend(loc='best', fontsize='small')
plt.axis([0,120,None,None])
plt.title('QMS current for H2, O2, Ar and Ar+O2 beams')
plt.ylabel('Scaled QMS current')
plt.xlabel('Time (min)')
plt.savefig('Figures/'+'allinone'+'.png',dpi=1000)
plt.show()
plt.close()