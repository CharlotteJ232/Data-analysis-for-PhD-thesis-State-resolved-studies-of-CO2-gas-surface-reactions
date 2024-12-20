from dataclasses import dataclass
from matplotlib import pyplot as plt
import matplotlib as mpl
import numpy as np
from scipy.optimize import least_squares, curve_fit
from brukeropusreader import read_file
from matplotlib.ticker import (AutoMinorLocator, MultipleLocator)
import os
from matplotlib import cm
import colorcet as cc
import math
import re


# mpl.rcParams['font.size']=20
textsize = 16
mpl.rcParams['xtick.labelsize']=textsize
mpl.rcParams['ytick.labelsize']=textsize
mpl.rcParams['axes.labelsize']=textsize




datalist = []

# # 20221205 With laser on resonance
# datalist.append({'folder': '2022/12 Dec/221205/RAIRS/',
# 'filestart': '7CO2_12He_100K_laseron',
# 'min_filenum': 0,
# 'max_filenum': 20 #max 20
# })

# # 20221205 With laser off resonance, dirty crystal
# datalist.append({'folder': '2022/12 Dec/221205/RAIRS/',
# 'filestart': '7CO2_12He_100K_laseroff',
# 'min_filenum': 0,
# 'max_filenum': 19 #max 19
# })

# # 20221205 With laser on resonance, dirty crystal
# datalist.append({'folder': '2022/12 Dec/221205/RAIRS/',
# 'filestart': '7CO2_12He_100K_laseron2',
# 'min_filenum': 0,
# 'max_filenum': 21 #max 21
# })

# # # 20221208 With laser off resonance
# datalist.append({'folder': '2022/12 Dec/221208/RAIRS/',
# 'filestart': '7CO2_12He_100K_laseroff',
# 'min_filenum': 0,
# 'max_filenum': 38 #max 39
# })

# # # 20221212 With laser on resonance
# datalist.append({'folder': '2022/12 Dec/221212/RAIRS/',
# 'filestart': '7CO2_12He_100K_laseron',
# 'min_filenum': 0,
# 'max_filenum': 28 #max 28
# })

# # # 20221213 With laser off resonance
# datalist.append({'folder': '2022/12 Dec/221213/RAIRS/',
# 'filestart': '7CO2_12He_100K_laseroff',
# 'min_filenum': 0,
# 'max_filenum': 20 #max 20
# })

# # 20221216 With laser off resonance
# datalist.append({'folder': '2022/12 Dec/221216/RAIRS/',
# 'filestart': '1CO2_14He_100K_laseroff',
# 'min_filenum': 0,
# 'max_filenum': 23 #max 23
# })

# # 20221219 With laser on resonance
# datalist.append({'folder': '2022/12 Dec/221219/RAIRS/',
# 'filestart': '1CO2_14He_100K_laseron',
# 'min_filenum': 0,
# 'max_filenum': 26 #max 26
# })

# # 20221219 Test if the water peak is still there after TPD
# datalist.append({'folder': '2022/12 Dec/221219/RAIRS/',
# 'filestart': 'test_water',
# 'min_filenum': 0,
# 'max_filenum': 0 #max 0
# })

# # 20221222 Oxygen dosing and then CO dosing
# datalist.append({'folder': '2022/12 Dec/221222/RAIRS/',
# 'filestart': '75s_O_16min_CO_Cu111',
# 'min_filenum': 0,
# 'max_filenum': 0 #max 0
# })

# # 20221223 Oxygen dosing and then CO dosing
# datalist.append({'folder': '2022/12 Dec/221223/RAIRS/',
# 'filestart': '120s_O_16mplus_CO_Cu111',
# 'min_filenum': 36,
# 'max_filenum': 36 #max 38
# })

# # # 20230102 Oxygen dosing at low Ekin and then CO dosing
# # datalist.append({'folder': '2023/01 Jan/230102/RAIRS/',
# # 'filestart': '4H_O2_0mplus_CO_Cu111',
# # 'min_filenum': 0,
# # 'max_filenum': 9 #max 9
# # })

# # 20230105 CO dosing on clean crystal, 5E-10 mbar
# datalist.append({'folder': '2023/01 Jan/230105/RAIRS/',
# 'filestart': '0m_O2_0mplus_CO_Cu111',
# 'min_filenum': 0,
# 'max_filenum': 6 #max 6
# })

# # 20230216 CO2 until saturation
# datalist.append({'folder': '2023/02 Feb/230216/RAIRS/',
# 'filestart': '4.75CO2_9He_100K_laseroff',
# 'min_filenum': 0,
# 'max_filenum':55  #max 55
# })

# # 20230221 H2 and then CO2
# datalist.append({'folder': '2023/02 Feb/230221/RAIRS/',
# 'filestart': 'D2_250K_CO2_100K',
# 'min_filenum': 23, #start CO2 dosing 23
# 'max_filenum':68  #max 68
# })

# # 20230302 H2 and then CO background dosing
# datalist.append({'folder': '2023/03 Mar/230302/RAIRS/',
# 'filestart': 'D2_250K_CO_180K',
# 'min_filenum': 0, #Double CO pressure starting from 11
# 'max_filenum':23  #max 23
# })

# # 20230327 Test measurements of just background stuff
# datalist.append({'folder': '2023/03 Mar/230327/RAIRS/',
# 'filestart': 'test_measurement_pausing',
# 'min_filenum': 0, 
# 'max_filenum':8  #max 8
# })

# # 20230327 Test measurements of just background stuff
# datalist.append({'folder': '2023/03 Mar/230328/RAIRS/',
# 'filestart': '2h_D2_long_CO2_Cu111',
# 'min_filenum': 1, 
# 'max_filenum':10  #max 10
# })

# # 20230522 CO buildup measurements with high Ekin beam and laser off
# datalist.append({'folder': '2023/05 May/230522/RAIRS/',
# 'filestart': '0.5CO2_9He_777C_100K_laseroff',
# 'min_filenum': 0, 
# 'max_filenum':35  #max 10
# })

# # 20240108 CO buildup measurements with high CO beam 
# datalist.append({'folder': '2024/01 Jan/240108/RAIRS/',
# 'filestart': 'CO_beam_100K',
# 'min_filenum': 0, 
# 'max_filenum':9  #max 9
# })

# # 20240108 CO buildup measurements with high CO beam 
# datalist.append({'folder': '2024/01 Jan/240108/RAIRS/',
# 'filestart': 'CO_beam_100K_afterflash',
# 'min_filenum': 0, 
# 'max_filenum':4  #max 4
# })

# # 20240108 CO buildup measurements with high CO beam 
# datalist.append({'folder': '2024/01 Jan/240108/RAIRS/',
# 'filestart': 'CO_beam_250K_afterflash',
# 'min_filenum': 1, 
# 'max_filenum':9  #max 9
# })

# # 240111 CO buildup measurements with high CO2 beam 
# datalist.append({'folder': '2024/01 Jan/240111/RAIRS/',
# 'filestart': '5CO2_250K',
# 'min_filenum': 0, 
# 'max_filenum':16  #max 19
# })

# # 240116 measurements with high CO2 beam and D2 coverage
# datalist.append({'folder': '2024/01 Jan/240116/RAIRS/',
# 'filestart': '5CO2_250K_withD2',
# 'min_filenum': 0, 
# 'max_filenum':16  #max 16
# })


# # 2024017 CO buildup measurements with CO beam 
# datalist.append({'folder': '2024/01 Jan/240117/RAIRS/',
# 'filestart': 'CO_beam_250K',
# 'min_filenum': 1, 
# 'max_filenum':9  #max 9
# })

# # 2024019 CO buildup measurements with CO beam and O2
# datalist.append({'folder': '2024/01 Jan/240119/RAIRS/',
# 'filestart': 'CO_beam_250K_withO2',
# 'min_filenum': 1, 
# 'max_filenum':6  #max 6
# })

# # 2024019 CO buildup measurements with CO beam and O2
# datalist.append({'folder': '2024/01 Jan/240119/RAIRS/',
# 'filestart': 'CO_beam_250K_withO2_afterflash',
# 'min_filenum': 1, 
# 'max_filenum':6  #max 6
# })

# 20240125 CO2 beam with D2
datalist.append({'folder': '2024/01 Jan/240125/RAIRS/',
#'filestart': 'CO2_beam_250K_withD2',
'filestart':'long',
'min_filenum': 0, 
'max_filenum':1  #max 16
})

# 20240223 CO beam with D2
datalist.append({'folder': '2024/02 Feb/240223/RAIRS/',
#'filestart': 'CO',
'filestart':'long',
'min_filenum': 0, 
'max_filenum':1  #22 long measurement, 23 after tpd
})



folderstart = "C:/Users/jansenc3/surfdrive/DATA/"
folderstart = "C:/Users/Werk/Surfdrive/DATA/"

use_every_x = 1 #only uses one in every x measurements
colors=['red', 'blue', 'lightcoral','cyan', 'purple', 'yellow']
save_data=False
save_figures=True
dpi = 300

yshift=0.007

def main():
    for dataset in datalist:
        load_dataset(dataset)

        npeak = 2078
        width = 10
        w_background = 3
        calculate_peak_properties(dataset, npeak, fit_voigt=False, width=width, w_background=w_background, save=save_data)

        # plot_peaks_single_dataset(dataset, npeak)

        dataset['spectrum_list'][-1].plot_raw(include_background=True, axis=[2115, 2035, 0.0225, None]) #plots last one only now
        print(dataset['spectrum_list'][-1].wavenumbers)
        print(len(dataset['spectrum_list'][-1].absorbance))

        if save_data: save_peaks_single_dataset(dataset, npeak)

        plot_spectra(dataset, use_every_x=use_every_x, yshift=yshift, nmin=1800, nmax=2400)
        plot_spectra(dataset, use_every_x=use_every_x, yshift=yshift, ymin=0, ymax=0.005)

        # for spectrum in dataset['spectrum_list']:
        #     spectrum.voigt_fit(npeak, plot=False)


    # plot_all_peaks_together(npeak)
    # # plot_integrals(npeak, typ='voigt')
    # plot_integrals(npeak, typ='sum')
    # plot_peak_position(npeak, typ='raw', yrange=5)
    # plot_peak_position(npeak, typ='voigt', yrange=5)
    # plot_peak_height(npeak)

    for dataset in datalist:
        npeak = 2050
        width = 250
        w_background = 10
        calculate_peak_properties(dataset, npeak, fit_voigt=False, width=width, w_background=w_background, save=save_data)

        # plot_peaks_single_dataset(dataset, npeak)

        if save_data: save_peaks_single_dataset(dataset, npeak)

    # plot_all_peaks_together(npeak)
    # plot_integrals(npeak, typ='sum')
    # # plot_peak_position(npeak, typ='raw', yrange=5)
    # # plot_integrals(npeak, typ='voigt')
    # # plot_peak_position(npeak, typ='voigt', yrange=5)
    # # plot_peak_height(npeak)

    for dataset in datalist:
        npeak = 2200
        width = 50
        w_background = 10
        calculate_peak_properties(dataset, npeak, width=width, w_background=w_background, save=save_data, fit_voigt=False)
        plot_firstlast(dataset,npeak,save=save_figures)
        plot_firstlastdiff(dataset,npeak,save=save_figures)

    plot_all_peaks_together(npeak)
    

    for dataset in datalist:
        npeak = 2850
        width = 200
        w_background = 10
        calculate_peak_properties(dataset, npeak, width=width, w_background=w_background, save=save_data, fit_voigt=False)
        plot_firstlast(dataset,npeak,save=save_figures)
        plot_firstlastdiff(dataset,npeak,save=save_figures)
   
    plot_all_peaks_together(npeak)

    for dataset in datalist:
        npeak = 1300
        width = 100
        w_background = 10
        calculate_peak_properties(dataset, npeak, width=width, w_background=w_background, save=save_data, fit_voigt=False)
        plot_firstlast(dataset,npeak,save=save_figures)
        plot_firstlastdiff(dataset,npeak,save=save_figures)
   
    plot_all_peaks_together(npeak)

    # for dataset in datalist:
    #     npeak = 2013
    #     width = 40
    #     w_background = 10
    #     calculate_peak_properties(dataset, npeak, width=width, w_background=w_background, save=save_data, fit_voigt=False)
    # plot_all_peaks_together(npeak)

    for dataset in datalist:
        npeak = 2000
        width = 200
        w_background = 20
        calculate_peak_properties(dataset, npeak, width=width, w_background=w_background, save=save_data, fit_voigt=False)
        plot_firstlast(dataset,npeak,save=save_figures)
        plot_firstlastdiff(dataset,npeak,save=save_figures)
    plot_all_peaks_together(npeak)


    for dataset in datalist:
        npeak = 3100
        width = 1900
        w_background = 20
        calculate_peak_properties(dataset, npeak, width=width, w_background=w_background, save=save_data, fit_voigt=False)
        plot_firstlastdiff_wide(dataset,npeak,save=save_figures, yrange=[-0.004, 0.003])
    plot_all_peaks_together(npeak)


############ Functions ###################

def load_dataset(dataset):
    min_filenum = dataset['min_filenum']
    max_filenum = dataset['max_filenum']
    folder = dataset['folder']
    filestart = dataset['filestart']
    
    dataset['spectrum_list'] = []
    filenumlist = np.arange(min_filenum,max_filenum+1,use_every_x)
    for filenum in filenumlist:
        if filenum < 10:
            filenumstr = '.000'+str(filenum)
        else:
            filenumstr = '.00'+str(filenum)
        file = filestart+filenumstr

        spectrum = Spectrum(filename=file)
        spectrum.load_data(folderstart+folder)
        dataset['spectrum_list'].append(spectrum)

def calculate_peak_properties(dataset, npeak, width=20, w_background=10, save=False, fit_voigt=True):
        dataset['times']=[]
        dataset['integrals']=[]
        dataset['peakpos']=[]
        dataset['peakheight']=[]
        dataset['voigt_integrals']=[]
        dataset['voigt_pos']=[]

        for spectrum in dataset['spectrum_list']:
            spectrum.get_time_from_data()
            dataset['times'].append(spectrum.timestamp)
            
            spectrum.simple_integral(npeak=npeak, width=width, w_background=w_background)
            dataset['integrals'].append(spectrum.peakdictionary[npeak]['integral'])
            
            if fit_voigt:
                spectrum.voigt_fit(npeak, plot=False)
                dataset['voigt_integrals'].append(spectrum.peakdictionary[npeak]['voigt_int'])
                dataset['voigt_pos'].append(spectrum.peakdictionary[npeak]['voigt_pos'])
                
            spectrum.simple_peakposition(npeak=npeak, width=width, w_background=w_background)
            dataset['peakpos'].append(spectrum.peakdictionary[npeak]['nmax'])
            dataset['peakheight'].append(np.max(spectrum.peakdictionary[npeak]['peak']))


        #Convert integral lists to np arrays and plot
        dataset['times'] = np.array(dataset['times'])/60 #convert to minutes
        dataset['integrals'] = np.array(dataset['integrals'])
        dataset['voigt_integrals'] = np.array(dataset['voigt_integrals'])
        dataset['peakpos'] = np.array(dataset['peakpos'])
        dataset['peakheight'] = np.array(dataset['peakheight'])
        dataset['voigt_pos']=np.array(dataset['voigt_pos'])

        if save:
            savefolder = folderstart+dataset['folder']+'Analyzed_data/'
            if not os.path.exists(savefolder):
                os.makedirs(savefolder)    
            np.savetxt(savefolder+dataset['filestart']+'_times.txt', dataset['times'].transpose(), header='Timestamp in minutes (from start of the day)')
            np.savetxt(savefolder+dataset['filestart']+'_integrals_'+str(npeak)+'.txt', dataset['integrals'].transpose(), header='Integrals of the peaks at '+str(npeak))
            np.savetxt(savefolder+dataset['filestart']+'_voigt_integrals_'+str(npeak)+'.txt', dataset['voigt_integrals'].transpose(), header='Integrals of the voigt peaks at '+str(npeak))
            np.savetxt(savefolder+dataset['filestart']+'_peakpos_'+str(npeak)+'.txt', dataset['peakpos'].transpose(), header='Position of the peaks at '+str(npeak))
            np.savetxt(savefolder+dataset['filestart']+'_peakheight_'+str(npeak)+'.txt', dataset['peakheight'].transpose(), header='Height of the peaks at '+str(npeak))
            np.savetxt(savefolder+dataset['filestart']+'_voigt_pos_'+str(npeak)+'.txt', dataset['voigt_pos'].transpose(), header='Position of the peaks at '+str(npeak)+' from voigt fit')



def save_peaks_single_dataset(dataset,npeak):
    data_tuple = (dataset['spectrum_list'][0].peakdictionary[npeak]['waven'],)
    for spectrum in dataset['spectrum_list']:
        data_tuple += (spectrum.peakdictionary[npeak]['peak'],)

    # print(data_tuple)
    data_array = np.column_stack(data_tuple)
    savefolder = folderstart+dataset['folder']+'Analyzed_data/'
    if not os.path.exists(savefolder):
        os.makedirs(savefolder)    
    np.savetxt(savefolder+dataset['filestart']+'_all_peaks_'+str(npeak)+'.txt', data_array, header='First column: wavenumbers, other columns: all peaks')

def plot_spectra(dataset, nmin=1200, nmax=5000, ymin=-0.01, ymax=0.05, yshift=0, start_spectrum=0, use_every_x=1):
    time_0 = dataset['spectrum_list'][0].timestamp
    fig, ax = plt.subplots(figsize=(6,4))
    for index in np.arange(start_spectrum, len(dataset['spectrum_list']), use_every_x):
        spectrum = dataset['spectrum_list'][index]
        ax.plot(spectrum.wavenumbers, spectrum.absorbance+index*yshift, linewidth=1, label=str(int((spectrum.timestamp-time_0)/60+0.5))+' min')
    ax.axis([nmax, nmin, ymin, ymax+index*yshift])
    ax.set_title(dataset['filestart'])
    ax.set_xlabel('Wavenumber')
    ax.set_ylabel('Absorbance')
    ax.legend(loc='upper left')
    
    plt.show()
    plt.close()

def plot_peak_height(npeak, save=False):
    print('plot_peak_height')
    fig, ax = plt.subplots()
    fig.set_dpi(dpi)
    for dataset, color in zip(datalist, colors): 
        times = dataset['times']-dataset['times'][0]          
        ax.plot(times, dataset['peakheight'],
                    'o', label=dataset['filestart'], color=color)
            


    ax.axis([0,None,None,None])
    ax.legend(loc='best')
    ax.set_xlabel('Time (minutes)')
    ax.set_ylabel('Peak height')

    savefolder = folderstart+dataset['folder']+'Analyzed_data/'
    if not os.path.exists(savefolder):
        os.makedirs(savefolder)   

def plot_peak_position(npeak, typ='raw', yrange=5):
    print('plot_peak_position')
    fig, ax = plt.subplots()
    fig.set_dpi(dpi)
    for dataset, color in zip(datalist, colors):   
        if typ == 'raw':
            data = dataset['peakpos']
        elif typ == 'voigt': 
            data = dataset['voigt_pos']

        ax.plot((dataset['times']-dataset['times'][0]), data,
                    'o', label=dataset['filestart'], color=color)
    ax.axis([0, None,npeak-yrange,npeak+yrange])
    ax.legend(loc='lower right')
    ax.set_title(typ)
    ax.set_xlabel('Time (minutes)')
    ax.set_ylabel('Peak position')

    fig, ax = plt.subplots()
    fig.set_dpi(dpi)
    for dataset, color in zip(datalist, colors):   
        if typ == 'raw':
            data = dataset['peakpos']
        elif typ == 'voigt': 
            data = dataset['voigt_pos']

        ax.plot(dataset['integrals'], data,
                    'o', label=dataset['filestart'], color=color)
    ax.axis([0, None,npeak-yrange,npeak+yrange])
    ax.legend(loc='lower right')
    ax.set_title(typ)
    ax.set_xlabel('Peak integral')
    ax.set_ylabel('Peak position')

def plot_peak_progression(spectrum_list, npeak, color, label=None):
    color_list = np.linspace(0, 0.9, len(spectrum_list))
    for spectrum, color in zip(spectrum_list, color_list):
            if not label:
                label=spectrum.filename[-2:]
            plt.plot(spectrum.peakdictionary[npeak]['waven'], spectrum.peakdictionary[npeak]['peak'],
                            color=[color,color,color], label=label)


def plot_all_peaks_together(npeak):
    print('plot_all_peaks_together')
    fig_peaks, ax_peaks = plt.subplots() 
    fig_peaks.set_dpi(dpi)
    for dataset, color in zip(datalist, colors):
        plot_peak_progression(dataset['spectrum_list'], npeak, color)
        # ax_peaks.plot([],[], color=color, label=dataset['filestart'])
    waven = dataset['spectrum_list'][0].peakdictionary[npeak]['waven']
    ax_peaks.axis([np.max(waven), np.min(waven), None,None])
    ax_peaks.set_xlabel('Wavenumber')
    ax_peaks.set_ylabel('Absorbance')
    # ax_peaks.legend()

def plot_firstlast(dataset,npeak,save=False, yrange=None):
    """
    does the same as plot_all_peaks together
    """
    print('plot_firstlast')
    fig_peaks, ax_peaks = plt.subplots() 
    fig_peaks.set_dpi(dpi)
    for spectrum, color, label in zip(dataset['spectrum_list'], ['black', 'blue'], ['after dosing', 'after TPD']):
        ax_peaks.plot(spectrum.peakdictionary[npeak]['waven'],spectrum.peakdictionary[npeak]['peak'], color=color, label=label)
        # ax_peaks.plot([],[], color=color, label=dataset['filestart'])
    waven = dataset['spectrum_list'][0].peakdictionary[npeak]['waven']
    if yrange:
            ax_peaks.axis([np.max(waven), np.min(waven), yrange[0],yrange[1]])
    else:
        ax_peaks.axis([np.max(waven), np.min(waven), None,None])
    ax_peaks.legend(loc='upper right')
    ax_peaks.set_xlabel('Wavenumber')
    ax_peaks.set_ylabel('Absorbance')
    if save:
        savefolder = folderstart+dataset['folder']+'Figures/'
        if not os.path.exists(savefolder):
            os.makedirs(savefolder)    
        plt.savefig(savefolder+str(npeak)+'_firstlast.png', dpi=300, bbox_inches='tight')
        plt.savefig(savefolder+str(npeak)+'_firstlast.eps', dpi=300, bbox_inches='tight')
        plt.savefig(savefolder+str(npeak)+'_firstlast.pdf', dpi=300, bbox_inches='tight')
    # ax_peaks.legend()

def plot_firstlastdiff(dataset,npeak,save=False, yrange=None):
    """
    does the same as plot_all_peaks together
    """
    print('plot_firstlastdiff')
    fig_peaks, axes = plt.subplots(2,1,sharex=True) 
    fig_peaks.subplots_adjust(hspace=0)
    ax_peaks = axes[0]
    ax_diff = axes[1]
    fig_peaks.set_dpi(dpi)
    for spectrum, color, label in zip(dataset['spectrum_list'], ['black', 'blue'], ['after dosing', 'after TPD']):
        ax_peaks.plot(spectrum.peakdictionary[npeak]['waven'],spectrum.peakdictionary[npeak]['peak'], color=color, label=label)
        # ax_peaks.plot([],[], color=color, label=dataset['filestart'])
    waven = dataset['spectrum_list'][0].peakdictionary[npeak]['waven']
    waven1 = dataset['spectrum_list'][1].peakdictionary[npeak]['waven']
    data0 = dataset['spectrum_list'][0].peakdictionary[npeak]['peak']
    data1 = dataset['spectrum_list'][1].peakdictionary[npeak]['peak']
    data = data0-data1
    ax_diff.plot(waven, data, color='red', label='Difference (dosing-TPD)')
    ax_diff.legend(loc='upper right')

    #print(waven == waven1)
    if yrange:
            ax_peaks.axis([np.max(waven), np.min(waven), yrange[0],yrange[1]])
    else:
        ax_peaks.axis([np.max(waven), np.min(waven), None,None])
    ylims = ax_peaks.get_ylim()
    ax_diff.set_ylim(ylims)
    ax_peaks.legend(loc='upper right')
    ax_peaks.set_xlabel('Wavenumber')
    ax_peaks.set_ylabel('Absorbance')
    if save:
        savefolder = folderstart+dataset['folder']+'Figures/'
        if not os.path.exists(savefolder):
            os.makedirs(savefolder)    
        plt.savefig(savefolder+str(npeak)+'_firstlastdiff.png', dpi=300, bbox_inches='tight')
        plt.savefig(savefolder+str(npeak)+'_firstlastdiff.eps', dpi=300, bbox_inches='tight')
        plt.savefig(savefolder+str(npeak)+'_firstlastdiff.pdf', dpi=300, bbox_inches='tight')
    # ax_peaks.legend()

def plot_firstlastdiff_wide(dataset,npeak,save=False, yrange=None):
    """
    does the same as plot_all_peaks together
    """
    print('plot_firstlastdiff_wide')
    fig_peaks, axes = plt.subplots(2,1,sharex=True, figsize=(12,4)) 
    fig_peaks.subplots_adjust(hspace=0)
    ax_peaks = axes[0]
    ax_diff = axes[1]
    fig_peaks.set_dpi(dpi)
    for spectrum, color, label in zip(dataset['spectrum_list'], ['black', 'blue'], ['after dosing', 'after TPD']):
        ax_peaks.plot(spectrum.peakdictionary[npeak]['waven'],spectrum.peakdictionary[npeak]['peak'], color=color, label=label)
        # ax_peaks.plot([],[], color=color, label=dataset['filestart'])
    waven = dataset['spectrum_list'][0].peakdictionary[npeak]['waven']
    waven1 = dataset['spectrum_list'][1].peakdictionary[npeak]['waven']
    data0 = dataset['spectrum_list'][0].peakdictionary[npeak]['peak']
    data1 = dataset['spectrum_list'][1].peakdictionary[npeak]['peak']
    data = data0-data1
    ax_diff.plot(waven, data, color='red', label='Difference (dosing-TPD)')
    ax_diff.legend(loc='upper right')

    #print(waven == waven1)
    if yrange:
            ax_peaks.axis([np.max(waven), np.min(waven), yrange[0],yrange[1]])
    else:
        ax_peaks.axis([np.max(waven), np.min(waven), None,None])
    ylims = ax_peaks.get_ylim()
    ax_diff.set_ylim(ylims)
    ax_peaks.legend(loc='upper right')
    ax_peaks.set_xlabel('Wavenumber')
    ax_peaks.set_ylabel('Absorbance')
    if save:
        savefolder = folderstart+dataset['folder']+'Figures/'
        if not os.path.exists(savefolder):
            os.makedirs(savefolder)    
        plt.savefig(savefolder+str(npeak)+'_firstlastdiff.png', dpi=300, bbox_inches='tight')
        plt.savefig(savefolder+str(npeak)+'_firstlastdiff.eps', dpi=300, bbox_inches='tight')
        plt.savefig(savefolder+str(npeak)+'_firstlastdiff.pdf', dpi=300, bbox_inches='tight')
    # ax_peaks.legend()

def plot_firstlast_wide(dataset,npeak,save=False, yrange=None):
    """
    does the same as plot_all_peaks together
    """
    print('plot_firstlast_wide')
    fig_peaks, ax_peaks = plt.subplots(figsize=(12,4)) 
    fig_peaks.set_dpi(dpi)
    for spectrum, color, label in zip(dataset['spectrum_list'], ['black', 'blue'], ['after dosing', 'after TPD']):
        ax_peaks.plot(spectrum.peakdictionary[npeak]['waven'],spectrum.peakdictionary[npeak]['peak'],linewidth=0.5, color=color, label=label)
        # ax_peaks.plot([],[], color=color, label=dataset['filestart'])
    waven = dataset['spectrum_list'][0].peakdictionary[npeak]['waven']
    if yrange:
            ax_peaks.axis([np.max(waven), np.min(waven), yrange[0],yrange[1]])
    else:
        ax_peaks.axis([np.max(waven), np.min(waven), None,None])
    ax_peaks.legend(loc='upper right')
    ax_peaks.set_xlabel('Wavenumber')
    ax_peaks.set_ylabel('Absorbance')
    if save:
        savefolder = folderstart+dataset['folder']+'Figures/'
        if not os.path.exists(savefolder):
            os.makedirs(savefolder)    
        plt.savefig(savefolder+str(npeak)+'_firstlast.png', dpi=300, bbox_inches='tight')
        plt.savefig(savefolder+str(npeak)+'_firstlast.eps', dpi=300, bbox_inches='tight')
        plt.savefig(savefolder+str(npeak)+'_firstlast.pdf', dpi=300, bbox_inches='tight')
    # ax_peaks.legend()

def plot_peaks_single_dataset(dataset, npeak):
    print('plot_peaks_single_dataset')
    fig_peaks, ax_peaks = plt.subplots(figsize=(4,4)) 
    fig_peaks.set_dpi(dpi)
    plot_peak_progression(dataset['spectrum_list'], npeak, None)
    waven = dataset['spectrum_list'][0].peakdictionary[npeak]['waven']
    ax_peaks.ticklabel_format(axis='y',style='sci', scilimits=(1,2))
    ax_peaks.axis([np.max(waven), np.min(waven), 0,None])
    ax_peaks.set_xlabel('Wavenumber')
    ax_peaks.set_ylabel('Absorbance')
    ax_peaks.set_title(dataset['filestart'])

def plot_integrals(npeak, typ='voigt'):
    """
    Plots integrals as a function of time for all datasets
    """
    print('plot_integrals')
    fig_int, ax_int = plt.subplots()
    fig_int.set_dpi(dpi)
    for dataset, color in zip(datalist, colors): 
        if typ == 'voigt':
            integrals = dataset['voigt_integrals'] 
        elif typ == 'sum':
            integrals = dataset['integrals']
        ax_int.plot((dataset['times']-dataset['times'][0]), integrals,
                    'o', label=dataset['filestart'], color=color)
    ax_int.axis([0, None, 0, None])
    ax_int.set_title(typ)
    ax_int.legend(loc='best')
    ax_int.set_xlabel('Time (minutes)')
    ax_int.set_ylabel('Integrated peak')


@dataclass
class Spectrum:
    filename: str = None
    beam_molecules: dict = None #example: {'CO2':5, 'He':15}
    dosing_time: float = None #total dosing time of the molecule of interest
    waiting_time: float = None #total time since the start of the experiment, can be used to account for desorption losses etc.

    def load_data(self, folderpath):
        self.folderpath = folderpath
        self.all_data = read_file(folderpath+self.filename)
        #calculate frequency array for plotting
        self.nmin = self.all_data['Fourier Transformation (Rf)']['HFQ']
        self.nmax = self.all_data['Fourier Transformation (Rf)']['LFQ']
        n_datapoints = len(self.all_data['AB'])
        self.wavenumbers = np.linspace(self.nmax, self.nmin, n_datapoints)
        self.absorbance = self.all_data['AB']
        self.background = self.all_data['ScRf']
        self.raw = self.all_data['ScSm']
        self.description = self.all_data['Sample']['SNM']

    def get_time_from_data(self):
        # digits = 111 + len(self.filename)
        tab = 6
        digits = 111 + (math.ceil(len(self.filename)/tab)) * tab
        
        match=(re.search('/',self.all_data['History']))
        digits = match.start()-4 #searches for the / character in the date, then -4 sets the index at the start of the year number
        # print(digits)

        # digits = 139
        # print(len(self.filename))
        # print(self.all_data['History'])
        dt = self.all_data['History'][digits:digits+20]
        # print(dt)
        time = dt[-9:]
        # print(time[:2],time[3:5],time[6:])
        self.timestamp = 3600*int(time[:2])+60*int(time[3:5])+int(time[6:])
        # print (self.timestamp)

    def set_datamask(self, datamask):
        self.datamask = datamask
    
    def plot_raw(self, include_background=True, axis=None):
        print('plot_raw')
        fig, ax = plt.subplots()
        fig.set_dpi(dpi)
        print(len(self.wavenumbers), len(self.absorbance), len(self.raw), len(self.background))
        ax.plot(self.wavenumbers, self.raw, color='r', label='Raw spectrum')
        if include_background:
            ax.plot(self.wavenumbers, self.background[4:-2], color='b', label='Background') #not sure why background is larger than wavenumbers and the rest. this works well enough..
            ax.legend()
        ax.axis([5000, 1200, 0, None])
        if axis:
            ax.axis(axis)
        ax.set_xlabel('Wavenumber')
        ax.set_ylabel('Reflected light intensity')
        plt.show()
        plt.close()



    def calculate_ab_background(self, datamask=None, typ='lin'):
        print('ab_background')
        if not self.datamask:
            if datamask:
                self.set_datamask(datamask)
            else:
                print('Please provide datamask')
        if typ == 'lin':
            a, b = np.polyfit(self.wavenumbers[datamask], self.absorbance[self.datamask], 1)
            self.ab_background=a*self.wavenumbers + b

    def simple_peakposition(self, npeak=2075,  width=20, w_background=10):
        try: self.peakdictionary[npeak]
        except AttributeError: self.isolate_peak(npeak=npeak, width=width, w_background=w_background)
        peak = self.peakdictionary[npeak]['peak']
        waven = self.peakdictionary[npeak]['waven']
        self.peakdictionary[npeak]['nmax'] = waven[np.argmax(peak)]


    def isolate_peak(self, npeak=2075, width=20, w_background=10, save=False):
        #check if the dictionary already exists
        try: self.peakdictionary
        except AttributeError: self.peakdictionary = {}
        
        #make wavenumber window
        peak_low = self.wavenumbers > npeak - width
        peak_high = self.wavenumbers < npeak + width
        peak = self.absorbance[peak_low*peak_high]
        peak_waven = self.wavenumbers[peak_low*peak_high]

        # #subtract background, old version
        # background = np.average(peak[:w_background]+peak[-w_background:])
        # # print(background)
        # peak -= background/2 #/2 because I added two arrays in bg calculation
        # # plt.plot(peak_freq, peak)

        #subtract background
        low = np.average(peak[:w_background])
        high = np.average(peak[-w_background:])
        peak -= np.linspace(low, high, len(peak))


        #save peak data
        self.peakdictionary[npeak] = {}
        self.peakdictionary[npeak]['peak'] = peak
        self.peakdictionary[npeak]['waven'] = peak_waven


    def simple_integral(self, npeak=2075, width=20, w_background=10, name=None):
        """
        Calculates the integral by isolating the peak, subtracting a 
        single value for the background and then taking the 
        sum of the result

        npeak: wavenumbers
        width: wavenumbers
        w_background: datapoints
        """ 
        try: self.peakdictionary[npeak]
        except AttributeError: self.isolate_peak(npeak=npeak, width=width, w_background=w_background)
        except KeyError: self.isolate_peak(npeak=npeak, width=width, w_background=w_background)

        self.peakdictionary[npeak]['integral'] = np.sum(self.peakdictionary[npeak]['peak'])
        # integrals.append(np.sum(peak))


    def voigt_fit(self, npeak, plot=False):
        def voigt(waven, A, f, waven_0, gamma_0):
            C = ((waven - waven_0)/gamma_0)**2
            def gauss(C, A, gamma_0):
                fac1 = A/gamma_0*np.sqrt(4*np.log(2)/np.pi)
                fac2 = np.exp(-4*np.log(2)*C)
                return fac1*fac2

            def lorentz(C, A, gamma_0):
                upper = 2*A/np.pi/gamma_0
                lower = 1+4*C
                return upper/lower

            # def return_voigt(C, A, gamma_0, f):
            return f*lorentz(C, A, gamma_0)+(1-f)*gauss(C, A, gamma_0)

            # return return_voigt()

        par_guess = [1,0.5,npeak,1]
        par_guess = [] #from diyu
        par_guess = [0.01, 0.5, npeak, 5]
        popt = curve_fit(voigt, self.peakdictionary[npeak]['waven'], 
                        self.peakdictionary[npeak]['peak'], 
                        p0=par_guess, maxfev=10000)

        # print(popt)
        smooth_waven = np.arange(np.min(self.peakdictionary[npeak]['waven']), np.max(self.peakdictionary[npeak]['waven']), 0.001)
        smooth_voigt = voigt(smooth_waven, *popt[0])

        self.peakdictionary[npeak]['voigt'] = voigt(self.peakdictionary[npeak]['waven'], *popt[0])
        self.peakdictionary[npeak]['voigt_int'] = np.sum(self.peakdictionary[npeak]['voigt'])
        index = np.argmax(smooth_voigt)
        self.peakdictionary[npeak]['voigt_pos'] = smooth_waven[index]

        if plot:
            fig_v, ax_v = plt.subplots()
            fig_v.set_dpi(dpi)
            ax_v.plot(self.peakdictionary[npeak]['waven'], self.peakdictionary[npeak]['voigt'])
            ax_v.plot(self.peakdictionary[npeak]['waven'], self.peakdictionary[npeak]['peak'])






if __name__ == '__main__':
    main()


