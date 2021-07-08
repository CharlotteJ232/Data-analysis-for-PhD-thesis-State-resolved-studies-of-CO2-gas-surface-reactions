# -*- coding: utf-8 -*-
"""
Created on Thu Feb 21 14:15:02 2019

@author: Charlotte
"""

from __future__ import division
import numpy as np
from matplotlib import pyplot as plt
import scipy.optimize

file1 = 'h2_co2'
file2 = 'o2_co2'
file3 = 'ar_co2_co'
file5 = 'ar_o2_co2_co_cont'
file4 = 'ar_o2_co2_co_disc'
file6 = 'ar_o2_co2_co'

filename = file6
ext = '.asc'

time, data1, data2 = np.loadtxt('Data/'+filename+ext, usecols=(4,5,6), 
                                skiprows=16, unpack=True)                                                                               
data = data1


#-----------------------------------------------------------------------------#                
                               
fit = np.polyfit(time, np.log(data), 1, w=np.sqrt(data)) 
fitdata = time*fit[0] + fit[1]
fitdata = np.exp(fitdata)   

fit2 = scipy.optimize.curve_fit(lambda t,a,b,c: a*np.exp(b*t)+c,  time,  data, 
                                p0=(np.exp(fit[1]), fit[0],0)) 
fitdata2 =  fit2[0][0]*np.exp(fit2[0][1]*time)+fit2[0][2]

print 'Y = '+str(fit2[0][0])+'exp('+str(fit2[0][1])+'t) + '+str(fit2[0][2])
plt.plot(time/3600, data, time/3600, fitdata2)
plt.show()
plt.close()