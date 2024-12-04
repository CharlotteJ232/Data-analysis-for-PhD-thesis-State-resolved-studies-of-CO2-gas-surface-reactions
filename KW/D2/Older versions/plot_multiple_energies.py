# -*- coding: utf-8 -*-
"""
Created on Fri Feb 28 15:05:09 2020

@author: jansenc3
"""

import numpy as np
from matplotlib import pyplot as plt

folderstart = 'P:/DATA/'
folderstart = 'D:/Surfdrive/DATA/'
year = '2020/'
subfolder = '/KW/Images/stickingprob_vs_position.txt'

days = [ '02 Feb/200214','03 Mar/200303','03 Mar/200305', '03 Mar/200310','03 Mar/200311']
energies = [ 0.8, 2.9, 1.2,0.7,4.0]
colors = [ 'r','c', 'black','g','y']

slopeA = []
slopeB = []
offset = []

minim = 18.8
center = 20
maxim = 21.2
for day,energy,c in zip(days,energies,colors):
    pos, data, error = np.loadtxt(folderstart+year+day+subfolder, unpack=True, usecols=(0,1,2), skiprows=1)
    
    if np.average(pos)>100:
        pos = pos/10
        
    left1 = pos > minim
    left2 = pos < center 
    left = left1*left2
    right1 = pos > center 
    right2 = pos < maxim
    right = right1*right2
    slopeleft, offsetleft = np.polyfit(pos[left], data[left], deg=1)
    slopeB.append(slopeleft)
    sloperight, offsetright = np.polyfit(pos[right], data[right], deg=1)
    slopeA.append(sloperight)
    offset.append((offsetleft+slopeleft*center+offsetright+sloperight*center)/2)
    
    
    color = str(1-(energy / 6.1))
    color = c
    plt.scatter(pos, data,label=str(energy)+' kJ/mol',color=color)
    plt.errorbar(pos, data, yerr=error, fmt='none', capsize=5,color=color)
    plt.plot(np.arange(minim,center,0.1), offsetleft+slopeleft*np.arange(minim,center,0.1),color=color)
    plt.plot(np.arange(center,maxim,0.1), offsetright+sloperight*np.arange(center,maxim,0.1),color=color)
plt.ylabel('Sticking probability')
plt.xlabel('Position (mm)')
plt.legend()
plt.savefig(folderstart+'D2 sticking probability as a function of step density and kinetic energy/slopesfit.png',dpi=1000)
plt.show()
plt.close()

X = np.column_stack((energies, slopeA, slopeB))
np.savetxt(folderstart+'D2 sticking probability as a function of step density and kinetic energy/slopes.txt',X,header='Energy (kJ/mol), slope A (/mm), slope B (/mm)')

plt.scatter(energies, np.absolute(slopeA), label='A type', c='b')
plt.scatter(energies, np.absolute(slopeB), label='B type', c='r')
plt.axis([0,7,0,None])
plt.title('slope (s0/position on crystal)')
plt.xlabel('Kinetic energy (kJ/mol)')
plt.ylabel('Slope (/mm)')
plt.savefig(folderstart+'D2 sticking probability as a function of step density and kinetic energy/slopesplot.png',dpi=1000)
plt.show()
plt.close()

plt.scatter(energies, offset, label='A type', c='black')
plt.axis([0,7,0,None])
plt.title('sticking probability at (111)')
plt.xlabel('Kinetic energy (kJ/mol)')
plt.ylabel('s0')
plt.savefig(folderstart+'D2 sticking probability as a function of step density and kinetic energy/s0_111.png',dpi=1000)
plt.show()
plt.close()

plt.scatter(np.log(energies), np.log(np.absolute(slopeA)), label='A type', c='b')
plt.scatter(np.log(energies), np.log(np.absolute(slopeB)), label='B type', c='r')
plt.title('loglog')
plt.xlabel('log Kinetic energy (kJ/mol)')
plt.ylabel('log Slope (/mm)')
plt.savefig(folderstart+'D2 sticking probability as a function of step density and kinetic energy/slopesplotloglog.png',dpi=1000)
plt.show()
plt.close()

plt.scatter((energies), np.log(np.absolute(slopeA)), label='A type', c='b')
plt.scatter((energies), np.log(np.absolute(slopeB)), label='B type', c='r')
plt.title('log')
plt.xlabel('Kinetic energy (kJ/mol)')
plt.ylabel('log Slope (/mm)')
plt.savefig(folderstart+'D2 sticking probability as a function of step density and kinetic energy/slopesplotlog.png',dpi=1000)
plt.show()
plt.close()