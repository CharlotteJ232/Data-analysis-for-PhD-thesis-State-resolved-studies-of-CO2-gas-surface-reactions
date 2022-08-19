import numpy as np
from matplotlib import pyplot as plt
from matplotlib.ticker import (AutoMinorLocator, MultipleLocator)


folderstart = "C:/Users/jansenc3/surfdrive/DATA/"
folderstart = "C:/Users/Werk/Surfdrive/DATA/"
folder = '2022/07 Jul/220718/KW/'
file = 'KW01'
# file = 'KW01_pressuregauge'
ext = '.txt'

save=True

time, data = np.loadtxt(folderstart+folder+file+ext, usecols=(0,1), unpack=True, skiprows=3)
time -= time[0]

totaltime = time[-1]

#remove overload values
index = data < 1E30
time = time[index]
data = data[index]

px = 1/plt.rcParams['figure.dpi']
fig, ax = plt.subplots(figsize=(totaltime*px, 200*px))
ax.set_xlabel('time (s)')
ax.set_ylabel('')
ax.plot(time, data)
ax.axis([0, None, 0, None])
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.xaxis.set_major_locator(MultipleLocator(200))
ax.grid()
if save:
    plt.savefig(folderstart+folder+file+'.png', dpi=200)

plt.show()
plt.close()