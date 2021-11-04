import numpy as np
from matplotlib import pyplot as plt


folder = 'DATA/2021/11 Nov/211102/TPD/'
directory = 'C:/Users/jansenc3/surfdrive/'
# directory = 'C:/Users/Werk/surfdrive/'
filename = 'KW01'
ext = '.txt'


for filename in ['tpd02', 'tpd03']:
    path = directory+folder+filename+ext
    time, temp, data = np.loadtxt(path,unpack=True,skiprows=3)
    plt.plot(time, data, label=filename)
    # plt.plot(time, temp, label=filename+' T')
plt.axis([None,None,None,None])
plt.legend()
plt.show()
plt.close()



for filename in ['kw02']:
    path = directory+folder+filename+ext
    time, data = np.loadtxt(path,unpack=True,skiprows=3)
    plt.plot(time, data, label=filename)
plt.axis([None,None,None,None])
plt.legend()
plt.show()
plt.close()



