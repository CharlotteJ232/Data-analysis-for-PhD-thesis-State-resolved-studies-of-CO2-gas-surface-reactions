import numpy as np
from matplotlib import pyplot as plt

folderstart = 'P:/Surfdrive/DATA/'
folderstart = 'C:/Users/Werk/surfdrive/DATA/'
loadfolder = folderstart+'Laser/Simulations/CH4/'
refdatapath = folderstart+"Laser/Simulations/CH4/RAP_fluence_curves.txt"

f_z0_list = [(0.254, 0.106, 'g'), (0.96, 0.57, 'r'), (1, 1, 'b')]

refdata = np.loadtxt(refdatapath,skiprows=2, usecols=(0,1,3,5,7),unpack=False)
for i in range(4):
    plt.plot(refdata[:,0],refdata[:,i+1], color='gray', label='Helen\'s tests')

datadic = {}
for f_z0 in f_z0_list[:3]:
    datadic[f_z0] = {}
    datadic[f_z0]['power'], datadic[f_z0]['population'] = np.loadtxt(loadfolder+'f_'+str(f_z0[0])+'_z_'+str(f_z0[1])+'.txt',skiprows=1, unpack=True)
    plt.plot(datadic[f_z0]['power'], datadic[f_z0]['population'], label=str(f_z0), c=f_z0[2])
plt.axis([0,500,0,1])
plt.legend()
plt.xlabel('Laser power (mW)')
plt.ylabel('Excited population')
plt.show()
plt.close()