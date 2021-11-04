# -*- coding: utf-8 -*-
"""
Created on Tue May  7 15:31:03 2019

@author: jansenc3
"""

import numpy as np
from matplotlib import pyplot as plt
import os
import glob

folder = 'P:/DATA/2019/05 May/190502/TOF/'
folder = 'P:/DATA/2019/08 Aug/190822/rga/'
directory = 'C:/Users/jansenc3/surfdrive/'
directory = 'C:/Users/Werk/surfdrive/'
folder_list=[directory+'DATA/2021/07 Jul/210722/TOF/Images/', directory+'DATA/2021/07 Jul/210723/TOF/Images/',directory+'DATA/2021/07 Jul/210726/TOF/Images/']
filename = 'HighresRGAfullrange'
ext = '.txt'

def main():
    plot_deltaE(folder_list)
    # plot_deltav(folder_list)
    plot_v_multiple(folder_list, vmin=0,vmax=None, pmax=6200)

def plot_deltaE(folder_list):
    for fold in folder_list:
        for subfolder in glob.glob(fold+'*'):
            for subsubfolder in glob.glob(subfolder+'/*'):
                p_nozzle = np.loadtxt(subsubfolder+'/Pnozzle.txt')
                deltaE = np.loadtxt(subsubfolder+'/deltE.txt')
                plt.plot(p_nozzle, deltaE, 'o', c='black')
    plt.xlabel('Nozzle pressure (mbar)')
    plt.ylabel('DeltaE/E')
    plt.axis([0,None,None,None])
    plt.show()
    plt.close()

def plot_deltav(folder_list):
    for fold in folder_list:
        for subfolder in glob.glob(fold+'*'):
            for subsubfolder in glob.glob(subfolder+'/*'):
                p_nozzle = np.loadtxt(subsubfolder+'/Pnozzle.txt')
                deltav = np.loadtxt(subsubfolder+'/deltv.txt')
                plt.plot(p_nozzle, deltav, 'o', c='black')
    plt.xlabel('Nozzle pressure (mbar)')
    plt.ylabel('Deltav/v')
    plt.axis([0,None,None,None])
    plt.show()
    plt.close()

def plot_v_multiple(folder_list,pmax=5200, vmin=0,vmax=None):
    for fold in folder_list:
        for subfolder in glob.glob(fold+'*'):
            for subsubfolder in glob.glob(subfolder+'/*'):
                v, yv = np.loadtxt(subsubfolder+'/velocity.txt',unpack=True, skiprows=1)
                p_nozzle = np.loadtxt(subsubfolder+'/Pnozzle.txt')
                color = p_nozzle/pmax
                plt.plot(v,yv, label=str(int(p_nozzle)), c=str(color))
    plt.xlabel('Velocity (m/s)')
    plt.ylabel('Distribution')
    plt.axis([vmin,vmax,0,None])
    plt.legend()
    plt.show()
    plt.close()



def plot_single():
    time, signal = np.loadtxt(folder+filename+ext, unpack=True)

    plt.plot(time, signal)
    #plt.axis([128,130,None,None])
    plt.show()
    plt.close()

main()