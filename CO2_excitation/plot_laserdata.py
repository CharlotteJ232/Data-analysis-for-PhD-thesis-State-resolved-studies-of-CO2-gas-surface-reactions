# -*- coding: utf-8 -*-
"""
https://matplotlib.org/3.2.1/tutorials/introductory/customizing.html
https://stackoverflow.com/questions/26106552/matplotlib-style-library-not-updating-when-mplstyle-files-added-deleted
https://matplotlib.org/3.2.1/tutorials/introductory/customizing.html#customizing-with-matplotlibrc-files
@author: Charlotte
"""
import datetime
import numpy as np
from matplotlib import pyplot as plt
# plt.style.reload_library()
# plt.style.use('voorbeeld')

folder = "C:/Users/Werk/Surfdrive/DATA/Power and wavelength data/"
date = "2020-10-23"

wl, lockin = np.loadtxt(folder+date+".txt", usecols=(3,6),unpack=True)

timearray = np.genfromtxt(folder+date+'.txt',dtype='str')[:,1]
print(timearray[10])

datetimearray = np.empty(len(timearray),dtype=str)
timestamparray = np.zeros(len(timearray))
for i in range(len(timearray)):
    datetimearray[i] = datetime.time(int(timearray[i][:2]),int(timearray[i][3:5]),int(timearray[i][6:8]))
    timestamparray[i] = 3600*int(timearray[i][:2]) + 60*int(timearray[i][3:5]) + int(timearray[i][6:8])


print (timearray[10],timestamparray[-1],wl[0],lockin[0])


fig, (ax1,ax2) = plt.subplots(2,sharex=True)
fig.subplots_adjust(hspace=0)
ax1.plot(timestamparray,wl,'-',linewidth=1,label='Laser wavelength',c='orange')
ax2.plot(timestamparray,lockin,'-',linewidth=1, label='PED signal from CO2')
plt.axis([57500,59500,None,None])
ax1.set_ylim(4252.7,4252.75)
ax2.set_ylim(3,5.5)
ticks = ax2.get_xticks()
ticks = ticks-ticks[0]
ax2.set_xticklabels(ticks.astype(int))

ax1.set_ylabel('Wavelength (nm)')
ax2.set_ylabel('Lock-in signal')
ax2.set_xlabel('Time (s)')
ax1.legend()
ax2.legend()
plt.show()
plt.close()