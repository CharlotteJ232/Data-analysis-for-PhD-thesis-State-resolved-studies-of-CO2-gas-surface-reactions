import pickle
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.patches as mpatches
from curlyBrace import curlyBrace
import sys
import os




folderstart = 'P:/Surfdrive/'
folderstart = 'C:/Users/Werk/surfdrive/'

sys.path.append(folderstart+'Python/')

from KW.process_KW_OO_positions import Measurement

filenamestart='KW'
folder = folderstart+'DATA/2020/03 Mar/200311/KW/Processed_data/'
savefolder = 'C:/Users/Werk/surfdrive/Proefschrift/Pictures etc/Methods/'
positions = np.arange(179, 220, 2)
t_flag1_open = 9 #s, total time flag 1 is open
t_open_flag2 = 3 #s, time after opening of flag 1
t_dose = 3 #s, dose time of the beam 
t_background = 4 #amount of seconds left and right of the measurement that are kept in the datasets. this can be chosen
xmin = -6
xmax = 10

save=True

def main():
    with open(folder+'data.pickle', 'rb') as f:
        # The protocol version used is detected automatically, so we do not
        # have to specify it.
        measurements = pickle.load(f)


    for position in [positions[9]]:
        plot_set(measurements[position], save=save)
        # plot_analyzed_data(measurements[position], save=False)
        plot_explanation(measurements[position],save=save)


def plot_set(measurement, save=False):
    col_max = len(list(measurement.datadict.keys()))
    scale = 1E12
    timeshift = t_open_flag2 + t_background
    fig, (ax1, ax2, ax3) = plt.subplots(3, sharex=True, figsize=(6,9))
    fig.subplots_adjust(hspace=0)
    for name in measurement.datadict.keys():
        measurement_number = float(name.replace(filenamestart,''))         
        ax1.plot(measurement.timedict[name]-timeshift,measurement.datadict[name]*scale,color=str(measurement_number/col_max*0.7), linewidth=1)
    # plt.plot(0,0,label='First',color=str(0.7/col_max))
    ax1.plot(-10,-10, label='Individual data',color=str(0.7))
    ax1.legend(loc='center right')



    ax2.plot(list(measurement.timedict.values())[0]-timeshift, measurement.sum*scale,'-',label='Summed data')
    ax2.legend(loc='center right')
    for ax in (ax1, ax2):
        ax.axis([xmin,xmax,0,None])
        ax.set_ylabel('QMS current (pA)')
        ax.plot([0,0],[0,scale], '--', label='Flag 2 open', color='black')
    ax2.set_xlabel('Time (s)')

    ax3.plot(list(measurement.timedict.values())[0]-timeshift, measurement.normalizeddata, label='Normalized KW data')
    # plt.plot(list(measurement.timedict.values())[0], measurement.baseline/np.max(measurement.sum), label='Baseline')
    ax3.plot(list(measurement.timedict.values())[0]-timeshift, measurement.fit, ':',color='black', label = 'Sticking dip fit')
    ax3.plot([0,0],[0,scale], '--', color='black', label='Flag 2 open')

    ax3.scatter(0, measurement.stickingprob, label = 'Initial sticking probability', color='black')
    ax3.legend(loc='center right')
    ax3.axis([xmin, xmax,0,1.1])
    ax3.set_xlabel('Time (s)')
    ax3.set_ylabel('Sticking probability')
    # plt.title('Position = '+str(measurement.position)+ ', Sticking probability = '+str(np.round(measurement.stickingprob,3)))


    for ax, letter in zip((ax1, ax2, ax3),('a', 'b', 'c')):
        ax.tick_params(top=True, direction='in')  
        ax.tick_params(right=True, labelleft=True, direction='in')
        ax.text(0.05, 0.82, letter+')', transform=ax.transAxes, fontsize=12)


    # plt.title(str(measurement.position))
    if save:
        if not os.path.exists(savefolder):
            os.makedirs(savefolder)
        plt.savefig(savefolder+'KW_analysis.png', dpi=500)
        plt.savefig(savefolder+'KW_analysis.pdf', dpi=500)
    plt.show()
    plt.close()

def plot_explanation(measurement, save=False):
    timeshift = t_open_flag2 + t_background
    scale = 1E9

    fig, ax = plt.subplots(figsize=(5,5))
    ax.plot(list(measurement.timedict.values())[0]-timeshift, measurement.sum*scale,'-',label='Summed data')
    ax.text(-5, 0.2, '(1)')
    ax.text(-t_open_flag2, 0.87, '(2)')
    ax.text(0, 0.87, '(3)')
    ax.text(t_dose, 0.87, '(4)')
    ax.text(t_flag1_open-t_open_flag2, 0.87, '(5)')
    # ax.text(8, 0.2, '(5)')

    # ax.text(1, 0.5, 'Sticking')
    # arrow = mpatches.Arrow(2, 0.55, 0, 0.15)
    # ax.add_patch(arrow)
    # ax.text(0,0.65, '{',rotation='vertical', fontsize=60, fontweight='ultralight')
    curlyBrace(fig, ax, [3,0.7], [0,0.7], 0.15, bool_auto=True, str_text='\n Sticking', color='black', lw=2, int_line_num=1)



    ax.set_xlabel('Time (s)')
    ax.set_ylabel('QMS current (nA)')
    ax.axis([xmin, xmax,0,1])
    ax.tick_params(top=True, direction='in')  
    ax.tick_params(right=True, labelleft=True, direction='in')
    
    if save:
        if not os.path.exists(savefolder):
            os.makedirs(savefolder)
        plt.savefig(savefolder+'KW_explanation.png', dpi=500)
        plt.savefig(savefolder+'KW_explanation.pdf', dpi=500)

    plt.show()
    plt.close()


def plot_analyzed_data(measurement, save=False):
    timeshift = t_open_flag2 + t_background
    fig, ax3 = plt.subplots()
    # plt.plot(list(measurement.timedict.values())[0], measurement.sum/np.max(measurement.sum),label='Original summed data (a.u.)')
    ax3.plot(list(measurement.timedict.values())[0]-timeshift, measurement.normalizeddata, label='Normalized KW data')
    # plt.plot(list(measurement.timedict.values())[0], measurement.baseline/np.max(measurement.sum), label='Baseline')
    ax3.plot(list(measurement.timedict.values())[0]-timeshift, measurement.fit, ':',color='black', label = 'Sticking dip fit')
    ax3.plot([0,0],[0,1], '--', color='black', label='Flag 2 open')

    ax3.scatter(0, measurement.stickingprob, label = 'Initial sticking probability', color='black')
    ax3.legend(loc='center right')
    ax3.axis([xmin, xmax,0,1])
    ax3.set_xlabel('Time (s)')
    ax3.set_ylabel('Sticking probability')
    # plt.title('Position = '+str(measurement.position)+ ', Sticking probability = '+str(np.round(measurement.stickingprob,3)))
    if save:
        if not os.path.exists(savefolder):
            os.makedirs(savefolder)
        plt.savefig(savefolder+'analyzed_'+str(measurement.position)+'.png', dpi=500)   
        plt.savefig(savefolder+'analyzed_'+str(measurement.position)+'.pdf', dpi=500)     
    plt.show()
    plt.close()


main()