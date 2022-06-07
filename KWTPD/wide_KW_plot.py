import numpy as np
from matplotlib import pyplot as plt


folderstart = "C:/Users/jansenc3/surfdrive/DATA/"
# folderstart = "C:/Users/Werk/Surfdrive/DATA/"
folder = '2022/04 Apr/220405/KW/'
file = 'KW01'
# file = 'KW01_pressuregauge'
ext = '.txt'

time, data = np.loadtxt(folderstart+folder+file+ext, usecols=(0,1), unpack=True, skiprows=3)
time -= time[0]

totaltime = time[-1]

px = 1/plt.rcParams['figure.dpi']
fig, ax = plt.subplots(figsize=(totaltime*px, 200*px))
ax.set_xlabel('time (s)')
ax.set_ylabel('')
ax.grid()
ax.plot(time, data)
plt.savefig(folderstart+folder+file+'.png', dpi=200)

plt.show()
plt.close()