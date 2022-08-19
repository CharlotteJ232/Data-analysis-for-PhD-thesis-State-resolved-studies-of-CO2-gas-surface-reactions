from matplotlib import pyplot as plt
import numpy as np
from scipy.optimize import least_squares
from brukeropusreader import read_file
from matplotlib.ticker import (AutoMinorLocator, MultipleLocator)
import os
from matplotlib import cm
import colorcet as cc

folderstart = "C:/Users/jansenc3/surfdrive/DATA/"
folderstart = "C:/Users/Werk/Surfdrive/DATA/"
oldfolder = '2022/07 Jul/220704/RAIRS/converted_data/'



min_filenum = 0
use_every_x = 1 

# 20220617 CO leak valve dosing
folder = '2022/06 Jun/220617/RAIRS/'
filestart = 'CO.00'
min_filenum = 0
max_filenum = 4
time_between_measurements = 5 #minutes
use_every_x = 1 #only uses one in every x measurements


# #20220704 test CO2
# folder = '2022/07 Jul/220704/RAIRS/raw_data/'
# filestart = '00_background.00'
# max_filenum = 19
# # file_b = filestart+'00.dpt'
# # freq_b, signal_b = np.loadtxt(folderstart+folder+file_b, unpack=True, delimiter=',')
#     # plt.plot(freq_b, signal_b)
#     # plt.show()
#     # plt.close()

# #20220708 CO2 in He, unstable beam
# folder = '2022/07 Jul/220708/RAIRS/'
# filestart = '1CO2_25He.00'
# max_filenum = 19

# # #20220714 No beam
# # folder = '2022/07 Jul/220714/RAIRS/'
# # filestart = 'No_beam.00'
# # max_filenum = 17

# # # 20220714 CO2 in Ar, with flag 2 open(?)
# # folder = '2022/07 Jul/220714/RAIRS/'
# # filestart = '1CO2_5Ar.00'
# # max_filenum = 10

# # #20220714 CO2 in Ar, with flag 2 closed(?) without flashing after previous measurement
# # folder = '2022/07 Jul/220714/RAIRS/'
# # filestart = '1CO2_5Ar.00'
# # max_filenum = 17
# # min_filenum = 11

# # 20220718 CO2 in He
# folder = '2022/07 Jul/220718/RAIRS/'
# filestart = '1CO2_25He.00'
# min_filenum = 0
# max_filenum = 24
# time_between_measurements = 5 #minutes
# use_every_x = 4 #only uses one in every x measurements


# # 20220718 CO2 in Ar
# folder = '2022/07 Jul/220722/RAIRS/'
# filestart = '1CO2_5Ar.00'
# min_filenum = 0
# max_filenum = 46
# time_between_measurements = 5 #minutes
# use_every_x = 1 #only uses one in every x measurements

# # 20220725 CO2 pure
# folder = '2022/07 Jul/220725/RAIRS/'
# filestart = '17CO2.00'
# min_filenum = 0
# max_filenum = 20
# time_between_measurements = 5 #minutes
# use_every_x = 4 #only uses one in every x measurements

# 20220728 CO2 in He
folder = '2022/07 Jul/220728/RAIRS/'
filestart = '5CO2_15He.00'
min_filenum = 23
max_filenum = 28
time_between_measurements = 5 #minutes
use_every_x = 1 #only uses one in every x measurements

# 20220728 CO2 in He, after dosing
folder = '2022/07 Jul/220728/RAIRS/'
filestart = '5CO2_15He.00'
min_filenum = 30
max_filenum = 46
time_between_measurements = 5 #minutes
use_every_x = 4 #only uses one in every x measurements

fmin = 2000
fmax = 2200
# fmax = 4000
# fmin = 3500 #for fit and plot
# fmax = 3800
fmin = 1200
fmax = 5000
shift_factor = 0.005 #to shift the different spectra vertically
plot_lines = [(2343, '  v3'), (2078, 'CO'), (3708, 'v1+v3'), (3600, '2v2+v3'), (2283, '13CO2'), (2382, 'LO')]

save=True

def main():

    if save:
        if not os.path.exists(folderstart+folder+'Images/'):
            os.makedirs(folderstart+folder+'Images/')


    plot_overlapping(ymax=0.3, ymin=-0.1)
    plot_shifted()



def plot_shifted():
    maxmin = int((max_filenum-min_filenum)/use_every_x) #number of plotted datasets
    fig, ax = plt.subplots(figsize=(6,(maxmin+1)))

    ytext = (maxmin + 1.1)*shift_factor

    for lineinfo in plot_lines:
        if (lineinfo[0] > fmin and lineinfo[0] < fmax):
            plot_line(ax, *lineinfo, ytext)

    shift=0
    for filenum in np.arange(min_filenum,max_filenum+1,use_every_x):
        if filenum < 10:
            filenumstr = '0'+str(filenum)
        else:
            filenumstr = str(filenum)
        file = filestart+filenumstr
        all_data = read_file(folderstart+folder+file)
        #calculate frequency array for plotting
        wavenmax = all_data['Fourier Transformation (Rf)']['HFQ']
        wavenmin = all_data['Fourier Transformation (Rf)']['LFQ']
        n_datapoints = len(all_data['AB'])
        freq = np.linspace(wavenmin, wavenmax, n_datapoints)
    
        #shift dataset in the plot
        index = (freq>fmin) * (freq<fmax)
        ax.plot(freq, all_data['AB']+shift-np.average(all_data['AB'][index]), label=str(filenum*time_between_measurements)+' min', color=dose_to_color(filenum*time_between_measurements)) #, color=dose_to_color(filenum*time_between_measurements)
        shift += shift_factor

    ax.axis([fmax, fmin, -shift_factor, (maxmin+1)*shift_factor])
    ax.yaxis.set_major_locator(MultipleLocator((shift_factor)))
    time_ax = ax.twinx()
    time_ax.set_ylabel('Dosing time (min)')
    time_ax.axis([None,None,min_filenum*time_between_measurements-time_between_measurements*use_every_x, max_filenum*time_between_measurements+time_between_measurements*use_every_x])
    time_ax.yaxis.set_major_locator(MultipleLocator((time_between_measurements*use_every_x)))
    time_ax.tick_params(direction='in')

    ax.set_xlabel('Wavenumber')
    ax.set_ylabel('Absorbance')
    figure_look(ax)
    plt.subplots_adjust(left=0.15, right=0.9, bottom=0.15, top=0.9)

    if save:
        filename = filestart+str(min_filenum)+'-'+str(max_filenum)+ '_'+str(fmin)+'-'+str(fmax)+'use_every_'+str(use_every_x)+'.png'
        plt.savefig(folderstart+folder+'Images/'+filename, dpi=500)
    plt.show()
    plt.close()

def plot_overlapping(ymax=0.15, ymin=-0.1):
    fig, ax = plt.subplots()
    ytext = ymax*1.02
    for lineinfo in plot_lines:
        if (lineinfo[0] > fmin and lineinfo[0] < fmax):
            plot_line(ax, *lineinfo, ytext)


    for filenum in np.arange(min_filenum,max_filenum+1,use_every_x):
        if filenum < 10:
            filenumstr = '0'+str(filenum)
        else:
            filenumstr = str(filenum)
        file = filestart+filenumstr
        all_data = read_file(folderstart+folder+file)
        #calculate frequency array for plotting
        wavenmax = all_data['Fourier Transformation (Rf)']['HFQ']
        wavenmin = all_data['Fourier Transformation (Rf)']['LFQ']
        n_datapoints = len(all_data['AB'])
        freq = np.linspace(wavenmin, wavenmax, n_datapoints)
    
        index = (freq>fmin) * (freq<fmax)
        ax.plot(freq, all_data['AB']-np.average(all_data['AB'][index]), label=str(filenum*time_between_measurements)+' min', color=dose_to_color(filenum*time_between_measurements)) #, color=dose_to_color(filenum*time_between_measurements)

    ax.axis([fmax, fmin, ymin, ymax])
    ax.legend()

    ax.set_xlabel('Wavenumber')
    ax.set_ylabel('Absorbance')
    figure_look(ax)
    plt.subplots_adjust(left=0.15, right=0.9, bottom=0.15, top=0.9)

    if save:
        filename = filestart+str(min_filenum)+'-'+str(max_filenum)+ '_'+str(fmin)+'-'+str(fmax)+'use_every_'+str(use_every_x)+'_overlap.png'
        plt.savefig(folderstart+folder+'Images/'+filename, dpi=500)
    plt.show()
    plt.close()

def dose_to_color(time, max_time=256, typ='cet_glasbey_light'):
    val = time/max_time
    if val>1:
        val==1
    return cm.get_cmap(typ)(val)

def plot_line(ax, wavenumber, text, ytext):
    ax.plot([wavenumber,wavenumber],[-100,100], '--', lw=1, color='gray')
    ax.text(wavenumber, ytext, text, size=7, rotation=45)

def figure_look(ax):
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(axis='x', top=True, direction='in')  #top=True, 
    ax.tick_params(axis='y', left=True, direction='in') #, labelleft=True, right=True, 
    ax.tick_params(which='minor', top=True, direction='in')  
    ax.tick_params(which='minor', left=True, right=False, direction='in')




### Old things ###

def fit_background(signal_b, signal):
    x = least_squares(fitfunction, 1, args=(signal_b, signal))
    return x.x[0]

def fitfunction(x, background=None, dataset=None):
    A = x[0]
    return dataset-A*background

def plot_old():
    fig, ax = plt.subplots(figsize=(6,9))
    for filenum in np.arange(1,max_filenum+1,1):
        if filenum < 10:
            filenumstr = '0'+str(filenum)
        else:
            filenumstr = str(filenum)
        file = filestart+filenumstr
        
        freq, signal = np.loadtxt(folderstart+folder+file, unpack=True, delimiter=',')

        index_for_fit = (freq_b>fmin) * (freq_b<fmax)
        A = fit_background(signal_b[index_for_fit], signal[index_for_fit])

        shift = filenum*shift_factor
        ax.plot(freq, signal/(A*signal_b)+shift)

    ax.axis([fmin, fmax, None,None])
    ax.set_xlabel('Wavenumber')
    ax.set_ylabel('Absorbance')
    plt.show()
    plt.close()

main()
