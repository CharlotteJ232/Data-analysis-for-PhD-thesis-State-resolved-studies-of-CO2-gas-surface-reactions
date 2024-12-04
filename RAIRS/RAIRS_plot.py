from matplotlib import pyplot as plt
import numpy as np
from scipy.optimize import least_squares
from brukeropusreader import read_file
from matplotlib.ticker import (AutoMinorLocator, MultipleLocator)
import matplotlib as mpl
import os
from matplotlib import cm
import colorcet as cc
import re


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

# 20221202 Test water and CO adsorption, dirty crystal
folder = '2022/12 Dec/221202/RAIRS/'
filestart = 'Test_100K_3.00'
min_filenum = 0
max_filenum = 39
use_every_x = 1 #only uses one in every x measurements

# 20221205 With laser on resonance
folder = '2022/12 Dec/221205/RAIRS/'
filestart = '7CO2_12He_100K_laseron.00'
min_filenum = 0
max_filenum = 20
use_every_x = 1 #only uses one in every x measurements

# # 20221205 With laser off resonance
# folder = '2022/12 Dec/221205/RAIRS/'
# filestart = '7CO2_12He_100K_laseroff.00'
# min_filenum = 0
# max_filenum = 19
# use_every_x = 1 #only uses one in every x measurements

# # 20221205 With laser on resonance
# folder = '2022/12 Dec/221205/RAIRS/'
# filestart = '7CO2_12He_100K_laseron2.00'
# min_filenum = 0
# max_filenum = 15
# use_every_x = 1 #only uses one in every x measurements

# # 20221208 With laser off resonance
# folder = '2022/12 Dec/221208/RAIRS/'
# filestart = '7CO2_12He_100K_laseroff.00'
# min_filenum = 0
# max_filenum = 25
# use_every_x = 1 #only uses one in every x measurements

# 20230508 CO2 dosing at 80K, pure beam, long measurement, used in paper
folder =  '2023/05 May/230508/RAIRS/'
filestart  = 'CO2_dosing.00'
min_filenum = 0
max_filenum = 29 #max 46
use_every_x = max_filenum - min_filenum
shift_factor = 0.2

use_files = [0, 1, 16]
use_colors = ['black', 'green', 'blue']
# (3550, fmax, fmin, 2550, ymin=-0.01, ymax=0.40, start_time=start_time, use_files=use_files, use_colors=use_colors, scale_first=10, shift_step=0.1)

# # 20230508 CO2 dosing at 80K, pure beam, clean Cu(111), short measurements
# folder =  '2023/05 May/230512/RAIRS/'
# filestart  = 'Laseron_01.00'
# min_filenum = 0
# max_filenum = 18 #max 18
# use_every_x = 3
# shift_factor = 0.03
# use_files = [0, 6, 12, 18]
# use_colors = ['black', 'green', 'gray', 'purple']
# text='Clean'
# #(3550, fmax, fmin, 2550, ymin=-0.01, ymax=0.5, start_time=start_time, use_files=use_files, use_colors=use_colors, scale_first=1, shift_step=0.1, text=text)

# # 20230508 CO2 dosing at 80K, pure beam, clean Cu(111), short measurements
# folder =  '2023/05 May/230512/RAIRS/'
# filestart  = 'Laseroff_03.00'
# min_filenum = 0
# max_filenum = 18 #max 18
# use_every_x = 3
# shift_factor = 0.03
# use_files = [0, 6, 12, 18]
# use_colors = ['black', 'green', 'gray', 'purple']
# text='Slightly dirty'

# # 20230508 CO2 dosing at 80K, pure beam, slightly dirty Cu(111), short measurements
# folder =  '2023/05 May/230512/RAIRS/'
# filestart  = 'Laseroff_09.00'
# min_filenum = 0
# max_filenum = 18 #max 18
# use_every_x = 3
# shift_factor = 0.03
# use_files = [0, 6, 12, 18]
# use_colors = ['black', 'green', 'gray', 'purple']
# text='Dirty'

# # 20230602 dosing at 80K, pure beam, dirty Cu(111)
# folder =  '2023/06 Jun/230602/RAIRS/'
# filestart  = 'Laseroff_03.00'
# min_filenum = 1
# max_filenum = 35 #max 38
# use_every_x = max_filenum-min_filenum
# shift_factor = 0.01

# # 20230622 CO2 dosing at 80K, pure beam, dirty Cu(111), long measurement
# folder =  '2023/06 Jun/230622/RAIRS/'
# filestart  = 'Test_KW.00'
# min_filenum = 0
# max_filenum = 32 #max 32
# use_every_x = 4
# shift_factor = 0.01
# use_files = [0, 6, 12, 18]
# use_colors = ['black', 'green', 'gray', 'purple']
# text='Slightly dirty'

# # 20230511 CO2 dosing at 80K, pure beam, slightly dirty Cu(111), long measurement
# folder =  '2023/05 May/230511/RAIRS/'
# filestart  = 'CO2_dosing.00'
# min_filenum = 0
# max_filenum = 17 #max 17
# use_every_x = 2
# shift_factor = 0.03
# use_files = [0, 6, 12, 17]
# use_colors = ['black', 'green', 'gray', 'purple']
# text=''



fmin = 2250
fmax = 2500
# fmax = 3400
# fmax = 2550
# # fmin = 3500 #for fit and plot
fmin = 2050
fmax = 3750
# fmin = 1200
# fmax = 5000
# shift_factor = 0.005 #to shift the different spectra vertically
plot_lines = [(2343, r'$\nu_3$'+'\n'), (2078, 'CO'), (3708, r'$\nu_1+\nu_3$'), (3600, r'$2\nu_2+\nu_3$'), (2283, '$^{13}$CO$_2$'), (2382, 'LO \n'),(2328, '----- X')]
# plot_lines = []
save=False

def main():

    if save:
        if not os.path.exists(folderstart+folder+'Images/'):
            os.makedirs(folderstart+folder+'Images/')

    start_time = get_start_time()
    plot_background()
    print('Plot overlapping')
    plot_overlapping(ymax=0.12, ymin=-0.02, start_time=start_time)
    plot_overlapping_broken_axis(3550, fmax, fmin, 2550, ymin=-0.01, ymax=0.5, 
                                start_time=start_time, use_files=use_files, 
                                use_colors=use_colors, scale_first=1, shift_step=0.1, text=text) # (3550, fmax, fmin, 2550, ymin=-0.01, ymax=0.40, start_time=start_time, use_files=use_files, use_colors=use_colors, scale_first=10, shift_step=0.1)
    print('Plot shifted')
    plot_shifted(start_time=start_time, bottom_space=2, top_space=5)
    # quick_integral()

def quick_integral():
    integrals = []
    filenumlist = np.arange(min_filenum,max_filenum+1,1)
    times = time_between_measurements*filenumlist
    times -= times[0]
    for filenum in filenumlist:
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

        #calculate integral
        line_freq = 2075 #wavenumbers
        linewidth = 20 #waveneumbers
        n_for_averaging = 10 #data points 

        #isolate peak
        peak_low = freq > line_freq - linewidth
        peak_high = freq < line_freq + linewidth
        peak = all_data['AB'][peak_low*peak_high]
        peak_freq = freq[peak_low*peak_high]

        # print("number of data points in peak:", len(peak))

        #subtract background
        background = np.average(peak[:n_for_averaging]+peak[-n_for_averaging:])
        # print(background)
        peak -= background/2 #/2 because I added two arrays in bg calculation
        plt.plot(peak_freq, peak)

        #calculate integral
        integrals.append(np.sum(peak))

    plt.xlabel('Wavenumber')
    plt.ylabel('Absorbance')
    plt.title(filestart)
    plt.show()
    plt.close()

    plt.plot(times, integrals, 'o')
    plt.xlabel('Time (minutes)')
    plt.ylabel('Peak integral')
    plt.title(filestart)
    plt.show()
    plt.close()




    

def get_time_from_data(data, unit='min'):
    """
    unit can be 'min', 'sec', 'hr' 
    """
    # digits = 111 + len(self.filename)
    # tab = 6
    # digits = 111 + (math.ceil(len(filename)/tab)) * tab
    
    match=(re.search('/',data['History']))
    digits = match.start()-4 #searches for the / character in the date, then -4 sets the index at the start of the year number
    # print(digits)

    # digits = 139
    # print(len(self.filename))
    # print(self.all_data['History'])
    dt = data['History'][digits:digits+20]
    # print(dt)
    time = dt[-9:]
    # print(time[:2],time[3:5],time[6:])
    timestamp = 3600*int(time[:2])+60*int(time[3:5])+int(time[6:])
    if unit == 'min':
        timestamp /= 60
    elif unit =='sec':
        timestamp = timestamp
    elif unit=='hr':
        timestamp /= 3600
    else:
        print('Please choose hours, minutes or seconds')
    # print (self.timestamp)
    return timestamp

def get_start_time(unit='min'):
    file = filestart+'00'
    all_data = read_file(folderstart+folder+file)
    return (get_time_from_data(all_data, unit=unit))


def plot_background():
    print('not finished')
    file = filestart+'14'
    all_data = read_file(folderstart+folder+file)
    print(all_data.keys())
    digits = 111 + len(file)
    print(all_data['History'][digits:digits+20])

    wavenmax = all_data['Fourier Transformation (Rf)']['HFQ']
    wavenmin = all_data['Fourier Transformation (Rf)']['LFQ']
    n_datapoints = len(all_data['AB'])
    freq = np.linspace(wavenmin, wavenmax, n_datapoints)
    plt.plot(freq, all_data['ScSm'])
    plt.axis([fmax, fmin, 0.008, 0.016])
    plt.show()
    plt.close()

def plot_shifted(start_time=0, bottom_space=1, top_space=5):
    maxmin = int((max_filenum-min_filenum)/use_every_x) #number of plotted datasets
    fig, ax = plt.subplots(figsize=(6,(maxmin+1)))

    shift=0
    first=True
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
        if first:
            firsttime = get_time_from_data(all_data) - start_time
            first=False
        shift = (get_time_from_data(all_data, unit='min') - start_time) * shift_factor
        index = (freq>fmin) * (freq<fmax)
        ax.plot(freq, all_data['AB']+shift-np.average(all_data['AB'][index]), label=str(filenum*time_between_measurements)+' min', color=dose_to_color(filenum*time_between_measurements)) #, color=dose_to_color(filenum*time_between_measurements)


    ytext = (maxmin + 1.1)*shift_factor
    ytext = (shift+top_space*shift_factor) * 1.01
    for lineinfo in plot_lines:
        if (lineinfo[0] > fmin and lineinfo[0] < fmax):
            plot_line(ax, *lineinfo, ytext)

    ax.axis([fmax, fmin, (firsttime-bottom_space)*shift_factor, shift+top_space*shift_factor])
    # ax.yaxis.set_major_locator(MultipleLocator((shift_factor)))
    time_ax = ax.twinx()
    time_ax.set_ylabel('Dosing time (min)')
    time_ax.axis([None, None, firsttime-bottom_space, shift/shift_factor+top_space])
    # time_ax.axis([None,None,min_filenum*time_between_measurements-time_between_measurements*use_every_x, max_filenum*time_between_measurements+time_between_measurements*use_every_x])
    # time_ax.yaxis.set_major_locator(MultipleLocator((time_between_measurements*use_every_x)))
    time_ax.tick_params(direction='in')

    ax.set_xlabel('Wavenumber')
    ax.set_ylabel('Absorbance')
    figure_look(ax)
    plt.subplots_adjust(left=0.15, right=0.9, bottom=0.15, top=0.9)

    if save:
        filename = filestart+str(min_filenum)+'-'+str(max_filenum)+ '_'+str(fmin)+'-'+str(fmax)+'use_every_'+str(use_every_x)
        plt.savefig(folderstart+folder+'Images/'+filename+'.png', dpi=500, bbox_inches='tight')
        plt.savefig(folderstart+folder+'Images/'+filename+'.eps', dpi=500, bbox_inches='tight')
    plt.show()
    plt.close()

def plot_overlapping_broken_axis(xmin1, xmax1, xmin2, xmax2, ymin=-0.1, ymax=0.15, start_time=0, use_files=None, use_colors=None, scale_first=1, shift_step=0, text=''):
    fig, axes = plt.subplots(ncols=2, gridspec_kw={'width_ratios': [1, 2.5]})
    fig.subplots_adjust(wspace=0.1)

    ytext = ymax*1.02

    for ax in axes:
        first=True
        shift = 0
        if not use_files:
            use_files = np.arange(min_filenum,max_filenum+1,use_every_x)
        for filenum, color_index in zip(use_files, range(len(use_files))):
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


            time = int(get_time_from_data(all_data, unit='min') - start_time + 0.5)
        
            index = (freq>fmin) * (freq<fmax)
            if use_colors:
                color = use_colors[color_index]
            else:
                color = dose_to_color(filenum*time_between_measurements)
            if first:
                factor = scale_first
                if scale_first!= 1:
                    label_extension = ' (x'+str(scale_first)+')'
                else:
                    label_extension = '' 
                first = False
            else:
                factor = 1
                

            ax.plot(freq, (all_data['AB']-np.average(all_data['AB'][index]))*factor+shift, label=str(time)+' min'+label_extension, color=color) #, color=dose_to_color(filenum*time_between_measurements)
            shift += shift_step

    axes[0].axis([xmax1, xmin1, ymin, ymax])
    axes[1].axis([xmax2, xmin2, ymin, ymax])

    #Create axis break markers
    d = 2  # proportion of vertical to horizontal extent of the slanted line
    kwargs = dict(marker=[(-1, -d), (1, d)], markersize=12,
                linestyle="none", color='k', mec='k', mew=1, clip_on=False)
    axes[0].plot([1, 1], [0, 1], transform=axes[0].transAxes, **kwargs)
    axes[1].plot([0, 0], [0, 1], transform=axes[1].transAxes, **kwargs)

    #Create labels and lines
    for lineinfo in plot_lines:
        if (lineinfo[0] > xmin1 and lineinfo[0] < xmax1):
            plot_line(axes[0], *lineinfo, ytext)
        if (lineinfo[0] > xmin2 and lineinfo[0] < xmax2):
            plot_line(axes[1], *lineinfo, ytext)

    axes[1].legend(fontsize=11)
    for ax in axes:
        figure_look(ax)
    axes[1].set_xlabel('Wavenumber')
    axes[0].spines['right'].set_visible(False)
    axes[1].spines['left'].set_visible(False)
    axes[1].tick_params(axis='y', left=False, labelleft=False, right=True) #, labelleft=True, right=True, 
    axes[1].tick_params(which='minor', left=False, right=True, direction='in')
    
    axes[0].set_ylabel('Absorbance')
    plt.subplots_adjust(left=0.15, right=0.9, bottom=0.15, top=0.9)
    axes[0].text(0.1, 0.9, text, transform=axes[0].transAxes, bbox=dict(boxstyle='square', facecolor='white', edgecolor='gray', alpha=0.7))

    if save:
        filename = filestart+str(min_filenum)+'-'+str(max_filenum)+ '_'+str(fmin)+'-'+str(fmax)+'use_every_'+str(use_every_x)+'_overlap_broken'
        plt.savefig(folderstart+folder+'Images/'+filename+'.png', dpi=500, bbox_inches='tight')
        plt.savefig(folderstart+folder+'Images/'+filename+'.eps', dpi=500, bbox_inches='tight')
    plt.show()
    plt.close()

def plot_overlapping(ymax=0.15, ymin=-0.1, start_time=0):
    fig, ax = plt.subplots()
    ytext = ymax*1.02
    for lineinfo in plot_lines:
        if (lineinfo[0] > fmin and lineinfo[0] < fmax):
            plot_line(ax, *lineinfo, ytext)


    first=True
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


        time = int(get_time_from_data(all_data, unit='min') - start_time + 0.5)
    
        index = (freq>fmin) * (freq<fmax)
        ax.plot(freq, all_data['AB']-np.average(all_data['AB'][index]), label=str(time)+' min', color=dose_to_color(filenum*time_between_measurements)) #, color=dose_to_color(filenum*time_between_measurements)

    ax.axis([fmax, fmin, ymin, ymax])
    ax.legend()

    ax.set_xlabel('Wavenumber')
    ax.set_ylabel('Absorbance')
    figure_look(ax)
    plt.subplots_adjust(left=0.15, right=0.9, bottom=0.15, top=0.9)

    if save:
        filename = filestart+str(min_filenum)+'-'+str(max_filenum)+ '_'+str(fmin)+'-'+str(fmax)+'use_every_'+str(use_every_x)+'_overlap'
        plt.savefig(folderstart+folder+'Images/'+filename+'.png', dpi=500, bbox_inches='tight')
        plt.savefig(folderstart+folder+'Images/'+filename+'.eps', dpi=500, bbox_inches='tight')
    plt.show()
    plt.close()

def dose_to_color(time, max_time=256, typ='cet_glasbey_light'):
    val = time/max_time
    if val>1:
        val==1
    return cm.get_cmap(typ)(val)

def plot_line(ax, wavenumber, text, ytext, fontsize=14):
    ax.plot([wavenumber,wavenumber],[-100,100], '--', lw=1, color='gray')
    ax.text(wavenumber, ytext, text, size=fontsize, rotation=45)

def figure_look(ax, fontsize=14):
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(axis='x', top=True, direction='in')  #top=True, 
    ax.tick_params(axis='y', left=True, direction='in') #, labelleft=True, right=True, 
    ax.tick_params(which='minor', top=True, direction='in')  
    ax.tick_params(which='minor', left=True, right=False, direction='in')
    mpl.rcParams.update({'font.size': fontsize})




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
