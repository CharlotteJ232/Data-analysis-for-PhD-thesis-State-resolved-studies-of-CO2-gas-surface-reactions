# -*- coding: utf-8 -*-
"""
Created on Mon Mar 18 15:55:09 2019

@author: jansenc3
"""

#reminder: data[

from datetime import datetime
import numpy as np
from matplotlib import pyplot as plt

filename = '0319/co2_co_he.txt'
#filename = 'test.txt'
start_analysis = datetime(2019,3,10,11,00,00)

start_analysis = datetime.timestamp(start_analysis)

start_time = np.loadtxt(filename, skiprows=1, max_rows=1)
names = np.loadtxt(filename,skiprows=2, max_rows=1, dtype=str, delimiter=',')

print (names[0])
data = np.loadtxt(filename, skiprows=3)

usedata = data[:,0]+start_time>start_analysis
#print (data[:,0])

newdata = data[usedata,:]
#print (newdata)


for k in [1,2,3]:
    #print (newdata[k,:])
    newdata[:,k] = newdata[:,k]-np.min(newdata[:,k])
    newdata[:,k] = newdata[:,k]/np.max(newdata[:,k])
    plt.scatter(newdata[:,0]/60,newdata[:,k],marker='.',s=0.02,label=names[k])
plt.xlabel('Time(min)')
plt.ylabel('Scaled QMS signal (a.u.)')
plt.axis([0,None,0,None])
plt.legend(loc='lower right')
plt.show()
plt.close()