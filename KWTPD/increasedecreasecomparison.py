# -*- coding: utf-8 -*-
"""
Created on Fri Aug 16 11:13:45 2019

@author: jansenc3
"""

import numpy as np
from matplotlib import pyplot as plt

folder = 'P:/DATA/2020/03 Mar/200311/KW/'
pos1, data1, error1 = np.loadtxt(folder+'Images/Increase/stickingprob_vs_position.txt',unpack=True, usecols=(0,1,2), skiprows=1)
pos2, data2, error2 = np.loadtxt(folder+'Images/Decrease/stickingprob_vs_position.txt',unpack=True, usecols=(0,1,2), skiprows=1)

plt.scatter(pos1, data1,label='Increase')
plt.errorbar(pos1, data1, yerr=error1, fmt='none', capsize=5)
plt.scatter(pos2, data2,label='Decrease')
plt.errorbar(pos2, data2, yerr=error2, fmt='none', capsize=5)
plt.legend()
plt.savefig(folder+'Images/updown.png',dpi=1000)
plt.show()
plt.close()

#def test(a, b, c, d, e):
#    print (a, b, c, d, e)
#   
#a = 0
#A=[0,1]
#B=[2,3]
#
#test(a,*A,*B)
#
#time3, data3, hd3 = np.loadtxt('P:/DATA/2020/02 Feb/200211/HD/HD03.txt',unpack=True, usecols=(0,1,2), skiprows=3)
#time4, data4, hd4 = np.loadtxt('P:/DATA/2020/02 Feb/200211/HD/HD04.txt',unpack=True, usecols=(0,1,2), skiprows=3)
#
#plt.plot(time3, data3, label='900K')
#plt.plot(time4,data4, label='1200K')
#plt.xlabel('Time')
#plt.legend()
#plt.title('D2 King&Wells')
#plt.axis([0,None,0,1.5E-10])
#plt.savefig('P:/DATA/2020/02 Feb/200211/HD/Images/KW.png', dpi=1000)
#plt.show()
#plt.close()
#
#plt.plot(time3, hd3, label='900K')
#plt.plot(time4,hd4, label='1200K')
#plt.xlabel('Time')
#plt.title('HD formation')
#plt.legend()
#plt.axis([0,None,0.011E-10,0.25E-10])
#plt.savefig('P:/DATA/2020/02 Feb/200211/HD/Images/HD.png', dpi=1000)
#plt.show()
#plt.close()

