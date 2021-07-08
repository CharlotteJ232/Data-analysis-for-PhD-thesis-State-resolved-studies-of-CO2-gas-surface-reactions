# -*- coding: utf-8 -*-
"""
Created on Tue May  5 13:33:59 2020

@author: Charlotte
"""

import numpy as np
from matplotlib import pyplot as plt

folderstart = 'P:/DATA/'
folderstart = 'D:/Surfdrive/DATA/'
folder = 'D2 sticking probability as a function of step density and kinetic energy/Models/Hard cube/'
energies = [3.38, 12.06, 16.88]
datadic = {}


for Ewell in energies:
    Es, datadic[Ewell] = np.loadtxt(folderstart+folder+'trapping_prob_'+str(Ewell)+'_kjm.txt',skiprows = 1, unpack=True)
    plt.plot(Es, datadic[Ewell], marker='o',label=str(Ewell))
plt.legend()
plt.axis([0,7,0,1])
plt.title('Trapping probability for different well depths (at 150K)')
plt.xlabel('Molecule kinetic energy (kJ/mol)')
plt.ylabel('Trapping probability')
plt.savefig(folderstart+folder+'trapping_prob_all.png')
plt.show()
plt.close()