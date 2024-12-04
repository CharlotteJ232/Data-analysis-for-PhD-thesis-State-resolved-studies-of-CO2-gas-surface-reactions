from matplotlib import pyplot as plt
import numpy as np
from scipy.optimize import least_squares
from brukeropusreader import read_file
from matplotlib.ticker import (AutoMinorLocator, MultipleLocator, MaxNLocator)
import os
from matplotlib import cm
import colorcet as cc


folderstart = "C:/Users/jansenc3/surfdrive/DATA/Pressures/"
# folderstart = "C:/Users/Werk/Surfdrive/DATA/Pressures/"

system = "Lionfish_HV"
day = "221120"
day = "221202"
day = "230323"
day = "220708"

gauge_names = {'Lionfish_HV': ['SSB1', 'SSB2', 'SSB3', 'Main', 'None', 'None']}
usegauges = (0,1,2,3)
gauge_names = {'Lionfish_HV': ['SSB1', 'SSB2', 'SSB3', 'None', 'None', 'Main']}
usegauges = (0,1,2,5)

xmin=0; xmax=None; ymin=None; ymax=None

hr = 3600 
xmin = 0*hr
xmax = 17*hr
# ymin = 8E-8
# ymax = 5E-7



file = folderstart+'Pressure_'+system+'/'+day+'_pressure.txt'

data = np.loadtxt(file, usecols=(2,3,4,5,6,7,8), unpack=True)

timestamp = data[0]
pressures = data[1:]
start_time = timestamp[0]
timestamp -= start_time
if xmax:
    window_min = timestamp > xmin
    window_max = timestamp < xmax
    window = window_min*window_max
else:
    window=np.full(len(timestamp), True)

fig, ax = plt.subplots()

for i in range(6):
    if np.any(np.array(list(gauge_names.keys())) == system):
        label = gauge_names[system][i]
    else:
        label = 'Gauge'+str(i+1)
    ax.plot(timestamp[window], pressures[i][window], label=label)


ax.axis([xmin, xmax, ymin, ymax])
fig.legend(loc='center right')
ax.set_ylabel('Pressure (mbar)')
ax.set_xlabel('Time (seconds)')
ax.set_yscale('log')
ax.set_title(system+' '+day)
plt.show()
plt.close()


fig_sep, axes = plt.subplots(nrows=len(usegauges), sharex=True)
fig_sep.subplots_adjust(hspace=0)

even=False
for j in range(len(usegauges)):
    i = usegauges[j]
    if np.any(np.array(list(gauge_names.keys())) == system):
        label = gauge_names[system][i]
    else:
        label = 'Gauge'+str(i+1)
    axes[j].plot(timestamp[window], pressures[i][window], label=label)
    # axes[j].set_yscale('lin')
    # axes[j].set_ylabel(label)
    axes[j].tick_params(which='major', right=True, top=True, direction='in')
    if even:
        axes[j].tick_params(which='major', labelleft=False, labelright=True)
    even = not even
    axes[j].tick_params(which='minor', left=False, labelleft=False, bottom=True, top=True)
    axes[j].yaxis.set_major_locator(MaxNLocator(3))
    axes[j].text(0.85, 0.8, label, transform=axes[j].transAxes)



axes[j].set_xlim(xmin, xmax)
axes[usegauges[0]].set_title(system+' '+day)
plt.show()
plt.close()