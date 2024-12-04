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
# from RAIRS_analysis_plot_functions import *
import math
import re


# mpl.rcParams['font.size']=20
textsize = 16
mpl.rcParams['xtick.labelsize']=textsize
mpl.rcParams['ytick.labelsize']=textsize
mpl.rcParams['axes.labelsize']=textsize


folderstart = "C:/Users/jansenc3/surfdrive/DATA/"
folderstart = "C:/Users/Werk/Surfdrive/DATA/"
datalist = []

min = 0
max = 18
max2 = 38
max2 = 18

#-------------- All start at the same dataset ----------------------


# # # 20230512 laser on and off CO2 buildup on cold crystal
# datalist.append({'folder': '2023/05 May/230512/RAIRS/',
# 'filestart': 'Laseron_02',
# 'min_filenum':min, 
# 'max_filenum':max  #max 18
# })

# datalist.append({'folder': '2023/05 May/230512/RAIRS/',
# 'filestart': 'Laseroff_03',
# 'min_filenum':min, 
# 'max_filenum':max  #max 18
# })

# datalist.append({'folder': '2023/05 May/230512/RAIRS/',
# 'filestart': 'Laseroff_04',
# 'min_filenum':min, 
# 'max_filenum':max  #max 18
# })

# datalist.append({'folder': '2023/05 May/230512/RAIRS/',
# 'filestart': 'Laseron_05',
# 'min_filenum':min, 
# 'max_filenum':max  #max 18
# })

# datalist.append({'folder': '2023/05 May/230512/RAIRS/',
# 'filestart': 'Laseron_06',
# 'min_filenum':min, 
# 'max_filenum':max  #max 18
# })

# datalist.append({'folder': '2023/05 May/230512/RAIRS/',
# 'filestart': 'Laseroff_07',
# 'min_filenum':min, 
# 'max_filenum':max  #max 18
# })

# datalist.append({'folder': '2023/05 May/230512/RAIRS/',
# 'filestart': 'Laseron_08',
# 'min_filenum':min, 
# 'max_filenum':max  #max 18
# })

# datalist.append({'folder': '2023/05 May/230512/RAIRS/',
# 'filestart': 'Laseroff_09',
# 'min_filenum':min, 
# 'max_filenum':max  #max 18
# })

datalist.append({'folder': '2023/06 Jun/230602/RAIRS/',
'filestart': 'Laseroff_03',
'min_filenum':min, 
'max_filenum':max2  #max 38
})

datalist.append({'folder': '2023/06 Jun/230602/RAIRS/',
'filestart': 'Laseron_04',
'min_filenum':min, 
'max_filenum':max2  #max 38
})

datalist.append({'folder': '2023/06 Jun/230602/RAIRS/',
'filestart': 'Laseron_05',
'min_filenum':min, 
'max_filenum':max2  #max 38
})

datalist.append({'folder': '2023/06 Jun/230602/RAIRS/',
'filestart': 'Laseroff_06',
'min_filenum':min, 
'max_filenum':max2  #max 38
})

use_every_x = 1 #only uses one in every x measurements
colors_varying = ['red', 'blue', 'lightcoral','cyan', 'purple', 'yellow', 'green', 'gray', 'black']
colors_rb = ['red', 'blue', 'blue', 'red', 'red', 'blue', 'red', 'blue', 'cyan', 'orange', 'orange', 'cyan']
colors_rb = ['red', 'blue', 'blue', 'red', 'red', 'blue', 'red', 'blue', 'blue', 'red', 'red', 'blue']
colors_rb = ['red', 'blue', 'lightblue', 'pink']
labels = ['laser on resonance', 'laser off resonance', 'laser off resonance', 'laser on resonance']
markershapes = ['o', 's', 's', 'o', 'o', 's', 'o', 's', 'v', 'p', 'p', 'v']
markershapes = ['o', 's', 's', 'o']
colors_saturation = np.linspace(0.5,1,len(colors_rb))
colors = [color_saturation*np.array(mpl.colors.to_rgb(color_rb)) for color_saturation, color_rb in zip(colors_saturation, colors_rb)]
colors = colors_rb
save_data=False
save_figures=False
dpi = 300
subfolder = 'first4/useevery5/'
use_every_x = 9
subfolder = 'first18useevery9/' 
# use_every_x = 1
# subfolder = ''

axis = [None, None, None, None]

# start_dataset = 2
# stop_dataset = 9 + 6 - 2

# datalist = datalist[start_dataset-2: stop_dataset-1]
# colors_rb = colors_rb[start_dataset-2: stop_dataset-1]
# markershapes = markershapes[start_dataset-2: stop_dataset-1]
# colors_saturation = colors_saturation[start_dataset-2: stop_dataset-1]


def main():

    for dataset in datalist:
        load_dataset(dataset)
        if save_figures:
            savefolder = folderstart+dataset['folder']+'Figures/'+subfolder
            if not os.path.exists(savefolder):
                os.makedirs(savefolder)    




    for dataset in datalist:
        npeak = 2345
        width = 5
        w_background = 2
        calculate_peak_properties(dataset, npeak, fit_voigt=False, width=width, w_background=w_background, save=save_data)

        # plot_peaks_single_dataset(dataset, npeak)

        # dataset['spectrum_list'][-1].plot_raw(include_background=True, axis=axis) #plots last one only now
        # print(dataset['spectrum_list'][-1].wavenumbers)
        # print(len(dataset['spectrum_list'][-1].absorbance))

        if save_data: save_peaks_single_dataset(dataset, npeak)

        # plot_spectra(dataset, use_every_x=use_every_x, yshift=0.005, nmin=1800, nmax=2400)
        # plot_spectra(dataset, use_every_x=use_every_x, yshift=0.005)

        # for spectrum in dataset['spectrum_list']:
        #     spectrum.voigt_fit(npeak, plot=False)



    plot_all_peaks_together(npeak, save=save_figures)
    # # plot_integrals(npeak, typ='voigt')
    plot_integrals(npeak, typ='sum', save=save_figures)
    # plot_peak_position(npeak, typ='raw', yrange=5)
    # plot_peak_position(npeak, typ='voigt', yrange=5)
    # plot_peak_height(npeak)


    # for dataset in datalist:
    #     npeak = 3709
    #     width = 7
    #     w_background = 3
    #     calculate_peak_properties(dataset, npeak, fit_voigt=False, width=width, w_background=w_background, save=save_data)

    #     # plot_peaks_single_dataset(dataset, npeak)

    #     if save_data: save_peaks_single_dataset(dataset, npeak)

    # plot_all_peaks_together(npeak, save=save_figures)
    # plot_integrals(npeak, typ='sum', save=save_figures, title=r'$\nu_1+\nu_3$')
    # # # plot_peak_position(npeak, typ='raw', yrange=5)
    # # # plot_integrals(npeak, typ='voigt')
    # # # plot_peak_position(npeak, typ='voigt', yrange=5)
    # # # plot_peak_height(npeak)

    # for dataset in datalist:
    #     npeak = 3601
    #     width = 7
    #     w_background = 3
    #     calculate_peak_properties(dataset, npeak, width=width, w_background=w_background, save=save_data, fit_voigt=False)
    # plot_all_peaks_together(npeak, save=save_figures)
    # plot_integrals(npeak, typ='sum', save=save_figures)

    # for dataset in datalist:
    #     npeak = 2284
    #     width = 5
    #     w_background = 2
    #     calculate_peak_properties(dataset, npeak, width=width, w_background=w_background, save=save_data, fit_voigt=False)
    # plot_all_peaks_together(npeak, save=save_figures)
    # plot_integrals(npeak, typ='sum', save=save_figures, title='$^{13}$CO$_2$')


    for dataset in datalist:
        npeak = 2380
        width = 60
        w_background = 10
        calculate_peak_properties(dataset, npeak, width=width, w_background=w_background, save=save_data, fit_voigt=False)
    plot_all_peaks_together(npeak, save=save_figures)
    plot_integrals(npeak, typ='sum', save=save_figures)

    # for dataset in datalist:
    #     npeak = 2013
    #     width = 40
    #     w_background = 10
    #     calculate_peak_properties(dataset, npeak, width=width, w_background=w_background, save=save_data, fit_voigt=False)
    # plot_all_peaks_together(npeak, save=save_figures)

    plot_spectra(dataset, nmin=2050, nmax=3750, ymin=-0.02, ymax=0.15, yshift=0, start_spectrum=16, use_every_x=10, save=save_figures)

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
            savefolder = folderstart+dataset['folder']+'Analyzed_data/'+subfolder
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

def plot_spectra(dataset, nmin=1200, nmax=5000, ymin=-0.01, ymax=0.05, yshift=0, start_spectrum=0, use_every_x=1, save=False):
    time_0 = dataset['spectrum_list'][0].timestamp
    fig, ax = plt.subplots(figsize=(6,4))
    for index in np.arange(start_spectrum, len(dataset['spectrum_list']), use_every_x):
        spectrum = dataset['spectrum_list'][index]
        ax.plot(spectrum.wavenumbers, spectrum.absorbance+index*yshift, linewidth=1, label=str(int((spectrum.timestamp-time_0)/60+0.5))+' min')
    ax.axis([nmax, nmin, ymin, ymax+len(dataset['spectrum_list'])*yshift])
    ax.set_title(dataset['filestart'])
    ax.set_xlabel('Wavenumber')
    ax.set_ylabel('Absorbance')
    ax.legend(loc='upper left')
    
    if save:
        savefolder = folderstart+dataset['folder']+'Figures/'+subfolder
        plt.savefig(savefolder+'spectra.png', dpi=dpi)
        plt.savefig(savefolder+'spectra.eps', dpi=dpi)

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
    """
    Plots all spectra of a single dataset, zoomed in on a single peak
    color as specified, lightness depends on number of the spectrum
    """
    lightness_list = np.linspace(0.5, 0.9, len(spectrum_list))
    color_list = [np.array(mpl.colors.to_rgb(color))*i for i in lightness_list]
    color_list = [color for i in lightness_list]

    # linestyles = ['-', '--', ':', '-.']
    # for spectrum, color, ls in zip(spectrum_list, color_list, range(len(linestyles))):
    #         plt.plot(spectrum.peakdictionary[npeak]['waven'], spectrum.peakdictionary[npeak]['peak'], ls=linestyles[ls],
    #                         color=color, label=spectrum.filename[-2:])
    first = True
    for spectrum, color in zip(spectrum_list, color_list):
            if not label:
                label = spectrum.filename[-2:]
            if first:
                plt.plot(spectrum.peakdictionary[npeak]['waven'], spectrum.peakdictionary[npeak]['peak'],
                            color=color, label=label, linewidth=1)
                first=False
            else:
                plt.plot(spectrum.peakdictionary[npeak]['waven'], spectrum.peakdictionary[npeak]['peak'],
                            color=color, linewidth=1)


def plot_all_peaks_together(npeak, save=False):
    """
    Plot all peaks of all datasets together, with colors varying between datasets 
    and lightness indicating the number of the spectrum in each dataset
    """
    print('plot_all_peaks_together')
    fig_peaks, ax_peaks = plt.subplots() 
    fig_peaks.set_dpi(dpi)
    for dataset, color, label in zip(datalist, colors, labels):
        plot_peak_progression(dataset['spectrum_list'], npeak, color, label)
        # ax_peaks.plot([],[], color=color, label=dataset['filestart'])
    waven = dataset['spectrum_list'][0].peakdictionary[npeak]['waven']
    ax_peaks.axis([np.max(waven), np.min(waven), None,None])
    ax_peaks.set_xlabel('Wavenumber')
    ax_peaks.set_ylabel('Absorbance')
    ax_peaks.legend(loc='best', fontsize=10)

    figure_look(ax_peaks)

    if save:
        savefolder = folderstart+dataset['folder']+'Figures/'+subfolder
        plt.savefig(savefolder+'_all_peaks_'+str(npeak)+'.png', bbox_inches='tight', dpi=dpi)
        plt.savefig(savefolder+'_all_peaks_'+str(npeak)+'.eps', bbox_inches='tight', dpi=dpi)




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

def plot_integrals(npeak, typ='voigt', title=None, save=False):
    """
    Plots integrals as a function of time for all datasets
    """
    print('plot_integrals')
    fig_int, ax_int = plt.subplots()
    fig_int.set_dpi(dpi)
    for dataset, color, marker, label in zip(datalist, colors, markershapes, labels): 
        if typ == 'voigt':
            integrals = dataset['voigt_integrals'] 
        elif typ == 'sum':
            integrals = dataset['integrals']
        ax_int.plot((dataset['times']-dataset['times'][0]), integrals,
                    'o', color=color, marker=marker, label=label)
    ax_int.axis([0, None, 0, None])
    if not title:
        ax_int.set_title(typ+' wavenumber '+str(npeak))
    else:
        ax_int.set_title(title)
    ax_int.legend(loc='best')
    ax_int.set_xlabel('Time (minutes)')
    ax_int.set_ylabel('Integrated peak')

    figure_look(ax_int)

    if save:
        savefolder = folderstart+dataset['folder']+'Figures/'+subfolder
        plt.savefig(savefolder+'_integrals_'+typ+'_'+str(npeak)+'.png', bbox_inches='tight', dpi=dpi)
        plt.savefig(savefolder+'_integrals_'+typ+'_'+str(npeak)+'.eps', bbox_inches='tight', dpi=dpi)


def figure_look(ax, fontsize=14):
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(axis='x', top=True, direction='in')  #top=True, 
    ax.tick_params(axis='y', left=True, right=True, direction='in') #, labelleft=True, right=True, 
    ax.tick_params(which='minor', top=True, direction='in')  
    ax.tick_params(which='minor', left=True, right=True, direction='in')
    mpl.rcParams.update({'font.size': fontsize})


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







main()


