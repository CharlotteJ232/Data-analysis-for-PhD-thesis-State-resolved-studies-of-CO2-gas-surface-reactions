# -*- coding: utf-8 -*-
"""
Created on Mon Mar 25 15:29:36 2019

@author: jansenc3
"""

from matplotlib import pyplot as plt 
import numpy as np

def plot_vs_temp(temp, data, label, title, save=False, filepath=None):
    plt.plot(temp, data, label=label)
    plt.title(title)
    plt.xlabel('Temperature (K)')
    plt.ylabel('Partial pressure (a.u.)')
    plt.legend()
    if save:
        plt.savefig(filepath,dpi=1000)
    plt.show()
    plt.close()   
    
def plot_integral_vs_beamtime(integralx,integrals,filepath):
    plt.scatter(integralx,integrals)
    plt.xlabel('Beam time (min)')
    plt.ylabel('Integral of peak (a.u.)')
    plt.title('Integral')
    plt.axis([0,None, 0, np.max(integrals)*1.1])
    plt.savefig(filepath,dpi=1000)
    plt.show()
    plt.close()
    
def plot_peak_edges(time, data, left, right, label, title=None, save=False, filepath=None):
    plt.plot(time, data, label=label)
    plt.scatter([time[left],time[right]],[data[left],data[right]],
                c='r',label='Peak edge definition')
    plt.title(title)
    plt.xlabel('Time (s)')
    plt.ylabel('Partial pressure (a.u.)')
    plt.legend()
    if save:
        plt.savefig(filepath,dpi=1000)
    plt.show()
    plt.close()  
