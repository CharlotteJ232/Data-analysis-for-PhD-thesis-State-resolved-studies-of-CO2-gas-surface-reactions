# -*- coding: utf-8 -*-
"""
Created on Tue May  7 15:31:03 2019

@author: jansenc3
"""

import numpy as np
from matplotlib import pyplot as plt

folder = 'P:/DATA/2019/05 May/190502/TOF/'
folder = 'P:/DATA/2019/08 Aug/190822/rga/'
filename = 'HighresRGAfullrange'
ext = '.txt'

time, signal = np.loadtxt(folder+filename+ext, unpack=True)

plt.plot(time, signal)
#plt.axis([128,130,None,None])
plt.show()
plt.close()