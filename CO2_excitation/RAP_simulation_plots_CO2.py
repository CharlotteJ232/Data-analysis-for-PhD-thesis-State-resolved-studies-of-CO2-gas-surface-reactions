import numpy as np
from matplotlib import pyplot as plt

folderstart = 'P:/Surfdrive/DATA/'
folderstart = 'C:/Users/Werk/surfdrive/DATA/'
loadfolder = folderstart+'Laser/Simulations/CO2/'

f_z0_list = [(0.2, 0.28, 'b', 'lens'), (0.2, 0, 'r', 'no lens')]

datadic = {}
for f_z0 in f_z0_list:
    datadic[f_z0] = {}
    datadic[f_z0]['power'], datadic[f_z0]['population'] = np.loadtxt(loadfolder+'f_'+str(f_z0[0])+'_z_'+str(f_z0[1])+'.txt',skiprows=1, unpack=True)
    

for f_z0_sublist in [f_z0_list]:
    for f_z0 in f_z0_sublist:   
        plt.plot(datadic[f_z0]['power'], datadic[f_z0]['population'], label=str(f_z0[3]), c=f_z0[2])
    plt.axis([0,50,0,1])
    plt.legend()
    plt.xlabel('Laser power (mW)')
    plt.ylabel('Excited population')
    plt.show()
    plt.close()
