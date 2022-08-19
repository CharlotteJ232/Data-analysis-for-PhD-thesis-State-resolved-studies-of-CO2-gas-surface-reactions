from matplotlib import pyplot as plt
import numpy as np
from scipy.optimize import least_squares
from brukeropusreader import read_file
from matplotlib.ticker import (AutoMinorLocator, MultipleLocator)
import os
from matplotlib import cm
import colorcet as cc


folderstart = "C:/Users/jansenc3/surfdrive/DATA/Pressures/"
# folderstart = "C:/Users/Werk/Surfdrive/DATA/Pressures"

system = "Lionfish_HV"
day = "220708"

gauge_names = {'Lionfish_HV': ['SSB1', 'SSB2', 'SSB3', 'Main', 'None', 'None']}

xmin=0; xmax=None; ymin=None; ymax=None

# xmin = 44000
xmax = 5000
# ymin = 1E-7
# ymax = 1E-6



file = folderstart+'Pressure_'+system+'/'+day+'_pressure.txt'

data = np.loadtxt(file, usecols=(2,3,4,5,6,7,8), unpack=True)

timestamp = data[0]
pressures = data[1:]
start_time = timestamp[0]



for i in range(6):
    if gauge_names[system]:
        label = gauge_names[system][i]
    else:
        label = 'Gauge'+str(i+1)
    plt.plot(timestamp-start_time, pressures[i], label=label)


plt.axis([xmin, xmax, ymin, ymax])
plt.legend(loc='best')
plt.ylabel('Pressure (mbar)')
plt.xlabel('Time (seconds)')
plt.yscale('log')
plt.show()
plt.close()