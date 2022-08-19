from xml.etree.ElementInclude import include
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.patches as mpatches
from curlyBrace import curlyBrace
import sys
import os
from matplotlib import cm
import matplotlib.ticker as ticker

folderstart = 'P:/Surfdrive/'
folderstart = 'C:/Users/Werk/surfdrive/'

filenamestart='KW'

datainfo = []
#datainfo.append([folder, scaling, ymax, xmin, xmax])
datainfo.append([folderstart+'DATA/2020/03 Mar/200310/TOF/', 0.02, 600, 0.65, 1.3])
datainfo.append([folderstart+'DATA/2021/07 Jul/210727/TOF/', 3.3, 5000, 0.15, 0.4])
datainfo.append([folderstart+'Data/2021/09 Sep/210930/TOF/', 3, 6000, 0.25, 0.5])
datainfo.append([folderstart+'DATA/2022/03 Mar/220321/TOF/', 1.1, 300, 0.35, 0.6])
datainfo.append([folderstart+'DATA/2022/04 Apr/220406/TOF25sccmHe1sccmO2/', 0.6, 2000, 0.35, 0.6])
datainfo.append([folderstart+'DATA/2022/04 Apr/220406/TOFb50sccmHe1sccmO2/', 1.2, 2000, 0.35, 0.6])
datainfo.append([folderstart+'DATA/2022/04 Apr/220406/TOFc40sccmHe1sccmO2/', 1, 2000, 0.35, 0.6])

save=False

extra_factor = 2

def main():
    for info in datainfo:
        print(info[0])
        velocity_slip_plot(*info)



def velocity_slip_plot(loadfolder, scaling, ymax, xmin, xmax):

    data_files = [exp_numbers_to_str(np.arange(7), 'TOF'), exp_numbers_to_str(np.arange(7,14), 'TOF')]

    n_plots = 7
    fig, axes = plt.subplots(n_plots, figsize=(6,12), sharex=True)
    fig.subplots_adjust(hspace=0)


    for i in range(n_plots):

        t, amp = np.loadtxt(loadfolder+data_files[1][i]+'.Asc',unpack=True)
        background_index = np.argwhere(t<xmin)[-20:]
        background = np.average(amp[background_index])
        axes[i].plot(t, extra_factor*scaling*(amp-background), color='blue', label=data_files[1][i], lw=1)          

        t, amp = np.loadtxt(loadfolder+data_files[0][i]+'.Asc',unpack=True)
        background_index = np.argwhere(t<xmin)[-20:]
        background = np.average(amp[background_index])
        axes[i].plot(t, amp-background, color='red', label=data_files[0][i], lw=1)


        axes[i].set_ylim(-100, ymax)
        axes[i].tick_params(labelleft=False)
        axes[i].legend()

        axes[i].tick_params(top=True, left=False, direction='in')

    axes[-1].set_xlim(xmin, xmax)
    axes[-1].set_xlabel('Time (ms)')
    fig.text(0.1, 0.5, 'QMS signal (arb. units)', ha='center', va='center', rotation='vertical')
    if save:
        plt.savefig(loadfolder+'Images/Velocity_slip.png', dpi=500)
    plt.show()
    plt.close()      


def exp_number_to_str(number, filenamestart=''):
       return filenamestart + str(int(number/10)) + str(int(number%10))
    
def exp_numbers_to_str(numbers, filenamestart=''):
    """numbers: array or list with numbers"""
    strs = []
    for i in range(len(numbers)):
        name = exp_number_to_str(numbers[i],filenamestart)
        strs.append(name)
    return strs

main()