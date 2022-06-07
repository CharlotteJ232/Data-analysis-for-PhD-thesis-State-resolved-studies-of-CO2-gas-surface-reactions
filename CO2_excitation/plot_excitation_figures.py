# -*- coding: utf-8 -*-
"""
https://matplotlib.org/3.2.1/tutorials/introductory/customizing.html
https://stackoverflow.com/questions/26106552/matplotlib-style-library-not-updating-when-mplstyle-files-added-deleted
https://matplotlib.org/3.2.1/tutorials/introductory/customizing.html#customizing-with-matplotlibrc-files
@author: Charlotte
"""
import datetime
import numpy as np
import os
from matplotlib import pyplot as plt
from matplotlib.ticker import (AutoMinorLocator, MultipleLocator)
from matplotlib import cm
from scipy.optimize import curve_fit
# plt.style.reload_library()
# plt.style.use('voorbeeld')

folderstart = "C:/Users/jansenc3/surfdrive/"
folderstart = "C:/Users/Werk/Surfdrive/"
folder = folderstart+"DATA/Power and wavelength data/"
simulationsfolder = folderstart+'DATA/Laser/Simulations/CO2/220517/'
simulationsfolder = folderstart+'DATA/Laser/Simulations/CO2/220603/'
savefolder = folderstart+'DATA/Rotational state distribution/Images/Newsimulations/Withrotations/'
savefolder = simulationsfolder
save=True
fun = 'dbl_exp'
cmap = 'viridis'

datainfo = { 'R0':{'year':'2022', 'month':'02', 'day':'24', 'min_power':0.005}, #was 0.0036
            'R2':{'year':'2022', 'month':'03', 'day':'04', 'min_power':0.0046},
            'R4':{'year':'2022', 'month':'03', 'day':'04', 'min_power':0.0041},
            'R6':{'year':'2022', 'month':'02', 'day':'25', 'min_power':0.005},
            'R8':{'year':'2022', 'month':'02', 'day':'25', 'min_power':0.005},
            'R10':{'year':'2022', 'month':'03', 'day':'03', 'min_power':0.005}
            }

transitions = ['R0', 'R2', 'R4', 'R6', 'R8', 'R10']
# transitions = ['R0', 'R2', 'R4']

months = {'01':'Jan', '02':'Feb', '03':'Mar', '04':'Apr', '05':'May', '06':'Jun', '07':'Jul', '08':'Aug', '09':'Sep', '10':'Oct', '11':'Nov', '12':'Dec'}



def main():
    if save:
        if not os.path.exists(savefolder):
            os.makedirs(savefolder)


    datadic = read_all()
    for transition in transitions:
        data = datadic[transition]
        datadic[transition] = normalize_shift_PED_signal(data, datainfo[transition]['min_power'], window_loss=0.1)
    # plot_single(datadic,'R2')
    # plot_all_normalized(datadic) 
    # plot_all_separately(datadic)
    # rotational_temperature(datadic)
    plot_all_data_simulations(datadic)






def read_all():
    datadic = {} 
    for transition in transitions:
        info = datainfo[transition]
        folder = (folderstart+'DATA/'+info['year']+'/'+info['month']+' '+months[info['month']]+'/'+info['year'][-2:]+info['month']+info['day']+'/Laser/Images/'+transition+'/')    
        datadic[transition] = {}
        datadic[transition]['power'], datadic[transition]['lockin'] = np.loadtxt(folder+'lockin_vs_power.txt', unpack=True)
        datadic[transition]['offresonance'] = np.loadtxt(folder+'lockin_offresonance.txt')
        datadic[transition]['popt'] = np.loadtxt(folder+'fitparameters_'+fun+'.txt')

        datadic[transition]['simulation_power'], datadic[transition]['simulation'] = np.loadtxt(simulationsfolder+transition+'_averaged.txt', unpack=True)
        # datadic[transition]['simulation_power'], datadic[transition]['simulation'] = np.loadtxt(simulationsfolder+transition+'/averaged.txt', unpack=True)

    return datadic

def plot_single(datadic, transition):
    data = datadic[transition]
    color = grayscale(transition, typ=cmap)
    fig, ax = plt.subplots()
    ax.plot(data['power']*1000, data['lockin'],'.',markersize=2, color=color, label='PED data '+transition)
    ax.plot(data['power'][data['sort']]*1000, dbl_exp(data['power'][data['sort']],*data['popt']), color='black', label='Fit')
    ax.plot([-100,100], [data['max'],data['max']], '--', color='black')
    ax.plot([-100,100], [data['offresonance'],data['offresonance']], '--', color='gray')
    ax.set_xlim(0, 35)
    ax.set_xlabel('Laser power (mW)')
    ax.set_ylabel('PED signal (V)')
    ax.legend(loc='center right')
    ax.tick_params(top=True, direction='in')  
    ax.tick_params(right=True, labelleft=True, direction='in')
    if save:
        plt.savefig(savefolder+transition+'.png', dpi=500)

def plot_all_normalized(datadic):    
    fig, ax = plt.subplots()
    for transition in list(datadic.keys()):
        color = grayscale(transition, typ=cmap)
        data = datadic[transition]
        ax.plot(data['power_shifted']*1000, data['lockin_normalized'],'.',markersize=2, color=color, label=transition)
    for transition in list(datadic.keys()):
        data = datadic[transition]
        ax.plot(data['power_shifted'][data['sort']]*1000, dbl_exp(data['power'][data['sort']],*data['popt'])/data['offresonance'],linewidth=1, color='black')
    ax.plot([-100,100], [1,1], '--', color='gray')
    ymax = 2.5
    ymin = 0.75
    ax.axis([-2, 40, ymin, ymax])
    ax.legend(loc='best')
    ax.set_xlabel('Laser power (mW)')
    ax.set_ylabel('PED signal (normalized)')

    ax.tick_params(top=True, direction='in')  
    ax.tick_params(right=True, labelleft=True, direction='in')
    ax.tick_params(which='minor', top=True, direction='in')  
    ax.tick_params(which='minor', right=True, direction='in')
    ax.yaxis.set_major_locator(MultipleLocator((ymax-ymin)/7))

    ax.text(0.025, 0.9, 'a)', transform=ax.transAxes)


    if save:
        plt.savefig(savefolder+'all_data_and_fits.png', dpi=500)
    plt.show()
    plt.close()

def plot_all_separately(datadic):
    for transition in transitions:
        data = datadic[transition]
        plt.plot(data['power_normalized']*1000, data['population'], '.',markersize=1, label=transition)
        plt.plot(data['simulation_power'], data['simulation'], color='black', label='Simulation')
        plt.xlabel('Laser power (mW)')
        plt.ylabel('Population')
        plt.legend()
        if save:
            plt.savefig(savefolder+transition+'_data_simulation.png', dpi=500)
        plt.show()

def plot_all_data_simulations(datadic):
    fig, (ax1, ax2, ax3) = plt.subplots(3,1,sharex='all',figsize=(5,9))
    fig.subplots_adjust(hspace=0)
    for transition in transitions:
        data = datadic[transition]
        color = grayscale(transition, typ=cmap)
        ax1.plot(data['power_normalized']*1000, data['population'], '.', color=color, markersize=1, label=transition)
        ax2.plot(data['power_normalized'][data['sort']]*1000, data['normalized_fit'][data['sort']], color=color, markersize=1, label=transition)
        ax3.plot(data['simulation_power'], data['simulation'], color=color, label=transition)

    labels=['PED data', 'Fit', 'Simulation']
    for ax, label,letter in zip([ax1, ax2, ax3], labels, 'abc'):
        ax.plot([-100,100],[1,1],'--',c='gray')
        ax.set_ylabel('Population')
        ax.set_ylim(0, 1.2)
        ax.tick_params(top=True, direction='in')  
        ax.tick_params(right=True, labelleft=True, direction='in')
        ax.text(0.025, 0.9, letter+')', transform=ax.transAxes)
        ax.text(0.4, 0.1, label, transform=ax.transAxes)

        # ax.xaxis.set_minor_locator(AutoMinorLocator())
        # ax.tick_params(which='minor', top=True, direction='in')  
        # ax.tick_params(which='minor', right=True, direction='in')
    ax2.tick_params(labelright=True, labelleft=False)
    ax2.yaxis.set_label_position("right")
    ax3.set_xlabel('Laser power (mW)')
    ax3.set_xlim(0, 40)
    ax3.legend(loc='lower right')






    if save:
        plt.savefig(savefolder+'comparison.png',dpi=500)

    plt.show()
    plt.close()

def rotational_temperature(datadic):
    J = []
    population = []
    colors = []
    for transition in transitions:
        J.append(int(transition[1:]))
        population.append((datadic[transition]['max']-datadic[transition]['offresonance'])/datadic[transition]['offresonance'])
        colors.append(grayscale(transition,typ=cmap))
    J = np.array(J)
    population = np.array(population)

    theta_r = 0.561 #for CO2
    T_guess = 5
    A_guess = 0.1
    vibrational_ground_state_fraction = 0.95

    def boltzmann(J, A, T):
        return A*(2*J+1)*np.exp(-1*(J*(J+1)*theta_r)/(T))

    popt = curve_fit(boltzmann, J,population,p0=[A_guess,T_guess])
    # print(popt[0])

    many_J = np.arange(0,1000,2)
    Z = np.sum(boltzmann(many_J, *popt[0]))/vibrational_ground_state_fraction

    fig, ax = plt.subplots()

    # plt.plot(J, population, 'o', label='Measured data')
    # plt.plot(J, boltzmann(J, A_guess, T_guess), label='Guess, T='+str(T_guess)+'K')
    smooth_J = np.arange(0,10, 0.1)
    ax.plot(smooth_J, boltzmann(smooth_J, *popt[0])/Z, '--', label='Fit, T='+str(np.round(popt[0][1],2))+'K', color='gray')
    for i in range(len(population)):
        ax.plot(J[i], population[i]/Z, 'o', markersize=8, color=colors[i])
    ax.plot([-1],[-1],'o',label='Measured population (normalized)',color='gray')
    plt.xlabel('J')
    plt.ylabel('Population')
    plt.legend(loc='lower right')

    ax.tick_params(top=True, direction='in')  
    ax.tick_params(right=True, labelleft=True, direction='in')
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(which='minor', top=True, direction='in')  
    ax.tick_params(which='minor', right=True, direction='in')

    # plt.title('Fitted T = '+str(np.round(popt[0][1],2))+' K')
    plt.axis([-0.5,10.5,0,0.25])


    ax.text(0.025, 0.9, 'b)', transform=ax.transAxes)

    if save:
        plt.savefig(savefolder+'rotational_temperature.png', dpi=500)
    plt.show()
    plt.close()



def grayscale(transition, typ='viridis'):
    shift = 4
    max = 16
    val = (int(transition[1:])+shift)/max
    if typ == 'gray':
        return str(val)
    elif typ == 'red':
        return (val, 0, 0)
    elif typ == 'green':
        return (0,val, 0)
    elif typ == 'blue':
        return (val/2,val/2,val)
    else:
        return cm.get_cmap(typ)(val)

def normalize_shift_PED_signal(data, min_power, window_loss = 0.1):
    """
    input: data dictionary, min_power corresponding to the dataset, window_loss is an estimate for the loss of light in the final window before the power meter
    Normalizes signal based on off-resonance signal. Shifts signal to correct for non-resonant laser light
    """
    data['power_shifted'] = data['power'] - min_power
    data['power_normalized'] = data['power_shifted'] / (1-window_loss)

    data['lockin_normalized'] = data['lockin'] / data['offresonance']
    data['population'] = (data['lockin'] - data['offresonance']) / (dbl_exp(1000000, *data['popt']) - data['offresonance'])

    data['normalized_fit'] = (dbl_exp(data['power'], *data['popt']) - data['offresonance']) / (dbl_exp(1000000, *data['popt']) - data['offresonance'])

    data['max'] = dbl_exp(1000000, *data['popt'])

    data['sort'] = np.argsort(data['power'])
    return data




def dbl_exp(power, A1, A2, b1, b2, c1, c2, d):
    return A1*(1-np.exp(b1*(power-c1)))+A2*(1-np.exp(b2*(power-c2)))+d

main()


### OLD ####
# folder = "C:/Users/Werk/Surfdrive/DATA/Power and wavelength data/"
# date = "2020-10-23"

# wl, lockin = np.loadtxt(folder+date+".txt", usecols=(3,6),unpack=True)

# timearray = np.genfromtxt(folder+date+'.txt',dtype='str')[:,1]
# print(timearray[10])

# datetimearray = np.empty(len(timearray),dtype=str)
# timestamparray = np.zeros(len(timearray))
# for i in range(len(timearray)):
#     datetimearray[i] = datetime.time(int(timearray[i][:2]),int(timearray[i][3:5]),int(timearray[i][6:8]))
#     timestamparray[i] = 3600*int(timearray[i][:2]) + 60*int(timearray[i][3:5]) + int(timearray[i][6:8])


# print (timearray[10],timestamparray[-1],wl[0],lockin[0])


# fig, (ax1,ax2) = plt.subplots(2,sharex=True)
# fig.subplots_adjust(hspace=0)
# ax1.plot(timestamparray,wl,'-',linewidth=1,label='Laser wavelength',c='orange')
# ax2.plot(timestamparray,lockin,'-',linewidth=1, label='PED signal from CO2')
# plt.axis([57500,59500,None,None])
# ax1.set_ylim(4252.7,4252.75)
# ax2.set_ylim(3,5.5)
# ticks = ax2.get_xticks()
# ticks = ticks-ticks[0]
# ax2.set_xticklabels(ticks.astype(int))

# ax1.set_ylabel('Wavelength (nm)')
# ax2.set_ylabel('Lock-in signal')
# ax2.set_xlabel('Time (s)')
# ax1.legend()
# ax2.legend()
# plt.show()
# plt.close()