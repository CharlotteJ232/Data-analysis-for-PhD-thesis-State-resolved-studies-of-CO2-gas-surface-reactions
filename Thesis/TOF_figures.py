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
folder = folderstart+'Data/2021/09 Sep/210930/TOF/'
savefolder = 'C:/Users/Werk/surfdrive/Proefschrift/Pictures etc/Setup_Methods/'

alpha = 50000 #determines width of example fitfunction for chopper convolution
fitparameters = [1,0,566.433566433567,1568.0914069668336,1E-10,151.4048824617757] #for the chopper plot


chopper_freq = 250
save=True
cmap = 'viridis'

def main():

    fig, axes = plot_raw_TOF(include_fit=False)
    for ax in axes:
        figure_style(ax)
    if save:
        plt.savefig(savefolder+'TOF_raw', dpi=500)

    fig, axes = plot_raw_TOF(include_fit=True)
    for ax in axes:
        figure_style(ax)
    if save:
        plt.savefig(savefolder+'TOF_fit', dpi=500)


    fig, ax = plot_chopper_function(include_TOF_shapes=False)
    figure_style(ax)
    figure_letter(ax, 'a')
    if save:
        plt.savefig(savefolder+'Chopper.png', dpi=500)

    fig, ax = plot_chopper_function(include_TOF_shapes=True, t_axis = 0.06)
    figure_style(ax)
    figure_letter(ax, 'b')
    if save:
        plt.savefig(savefolder+'Chopper_realistic.png', dpi=500)


    global fitparameters 
    fitparameters[-1] = alpha
    fig, ax = plot_chopper_function(include_TOF_shapes=True, t_axis=0.015)
    figure_style(ax)
    figure_letter(ax, 'c')
    if save:
        plt.savefig(savefolder+'Chopper_example.png', dpi=500)



def figure_letter(ax, letter):
    ax.text(0.05, 0.9, letter+')', transform=ax.transAxes, fontsize=12)



def plot_raw_TOF(include_fit=False):
    fig, axes = plt.subplots(2, sharex=True)
    fig.subplots_adjust(hspace=0)
    masses = ['4.0', '44.0']
    nrs = [0,7]
    names = ['He','CO$_2$']
    if include_fit:
        data_files=[[exp_number_to_str(nrs[0], 'TOF')], [exp_number_to_str(nrs[1], 'TOF')]]
    else:
        data_files = [exp_numbers_to_str(np.arange(7), 'TOF'), exp_numbers_to_str(np.arange(7,14), 'TOF')]

    for i in range(2):
        colorcounter = 0
        for data_file in data_files[i]:
            if include_fit:
                color = 'lightgray'
            else:
                color = cm.get_cmap(cmap)(colorcounter/len(data_files[i]))
                colorcounter += 1
            t, amp = np.loadtxt(folder+data_file+'.Asc',unpack=True)
            background = np.average(amp[:20])
            axes[i].plot(t, amp-background, color=color, label='Raw data')
            if include_fit:
                filename = 'fit_TOF'+str(nrs[i])+'.0zoom.txt'
                tfit, yfit = np.loadtxt(folder+'Images/'+data_file+'/'+masses[i]+'/'+filename, usecols=(0,1), skiprows=1, unpack=True)
                axes[i].plot(tfit, yfit-background, color='black', label='Fit')
                axes[i].legend(loc='upper left')

        axes[i].axis([0.25, 0.55, None, None])
        axes[i].tick_params(labelleft=False)
        axes[i].text(0.9, 0.8, names[i], transform=axes[i].transAxes, fontsize=12)
    
    fig.text(0.1, 0.5, 'QMS signal (arb. units)', ha='center', va='center', rotation='vertical')
    
    axes[1].set_xlim(0.25, 0.55)
    axes[0].set_ylim(-100, 6000)
    axes[1].set_ylim(-100, 1500)
    axes[1].set_xlabel('Time (ms)')

    return fig, axes

def exp_number_to_str(number, filenamestart=''):
       return filenamestart + str(int(number/10)) + str(int(number%10))
    
def exp_numbers_to_str(numbers, filenamestart=''):
    """numbers: array or list with numbers"""
    strs = []
    for i in range(len(numbers)):
        name = exp_number_to_str(numbers[i],filenamestart)
        strs.append(name)
    return strs

def figure_style(ax):
    ax.tick_params(top=True, right=True, direction='in')

def plot_chopper_function(include_TOF_shapes=False, t_axis=0.015):
    smooth_chopper = get_chopper_amplitudes(n=200, chopper_freq=chopper_freq)
    fig, ax = plt.subplots()
    ax.plot(smooth_chopper[0], smooth_chopper[1], linewidth=2, color='black', label='Beam pulse')
    ax.plot([0,0],[0,1], '--', color='gray', label='Delta function')
    ax.set_xlabel('Time ($\mu$s)')
    ax.set_ylabel('Beam amplitude (relative)')

    if include_TOF_shapes:
        plot_shifted_scaled_fitfunctions(ax, fitparameters, n=3, chopper_freq=chopper_freq)

    tscale = 1000
    ax.axis([-t_axis, t_axis, 0, 1.1])
    ticks = ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(x*tscale))
    ax.xaxis.set_major_formatter(ticks)
    ax.legend(loc='upper right', fontsize=9)

    return fig, ax


def plot_shifted_scaled_fitfunctions(ax, fitparameters, n, chopper_freq=250):
    tchop_list, amp_list = get_chopper_amplitudes(n=n, chopper_freq=chopper_freq)
    twindow = 2
    # tcenter = fitparameters[-4]/fitparameters[-3]
    tcenter=0
    t = np.arange(tcenter-twindow, tcenter+twindow,0.00001)

    ylist = fitfunction_convolution(t, *fitparameters, n=n, chopper_freq=chopper_freq)[0]
    tchop, achop = get_chopper_amplitudes(n=n, chopper_freq=chopper_freq)
    print(tchop)
    tcenter = t[np.argmax(ylist[int(len(ylist)/2)])]
    print(tcenter)
    for i in range(n):
        tcutoff = t > tchop[i]
        ax.plot([tchop[i], tchop[i]],[0,achop[i]], '--', color='gray')
        ax.plot(t[tcutoff]-tcenter, ylist[i][tcutoff]/np.max(ylist[i])*achop[i], color=cm.get_cmap(cmap)(i/n), label='TOF function')




def fitfunction_convolution(t, A, B, L, t0, ts, alpha, chopper_freq=250,corr_dens=False, n=15, pos=0, use_v0=True):
    """
    A = amplitude
    B = background
    L = flight path length
    t0 = peak maximum position
    ts = time shift of the function (=delay time)
    alpha = width of the peak
    
    chopper_freq = frequency of the chopper
    corr_dens = correct for density sensitive detector (to get the actual shape 
    of the peak in the time domain, instead of what is measured)
    n = number of points for convolution
    pos = shift in L, position of mass spectrometer
    """
    t_shift, amp = get_chopper_amplitudes(n, chopper_freq) # shift/amplitude array for convolution

    L = L-pos  
    
    if use_v0:
        t0 = L/t0 #because given t0 was actually v0
  
    y = np.zeros(len(t))
    ylist = []
    for i in range(n):
        #prevent dividing by 0
        # import pdb; pdb.set_trace()
        if np.any(t-ts-t_shift[i]==0):
            print ('dividing by 0')
            return 0

        if corr_dens:
            #note: i think the last /t should instead be / (t-ts-t_shift[i]), i don't know why I chose for just /t.
            #also, this version cannot be used to fit, only to generate the y values for certain fit parameters
            ysep = (amp[i] * A * ((L/(t-ts-t_shift[i]))**4 * np.exp(-((L/(t-ts-t_shift[i])-L/t0)/alpha)**2)) / t)
            ylist.append(ysep)
            y += ysep #
        else:
            ysep = (amp[i] * A * ((L/(t-ts-t_shift[i]))**4 * np.exp(-((L/(t-ts-t_shift[i])-L/t0)/alpha)**2))) #
            ylist.append(ysep)
            y += ysep
    y /= np.sum(amp)  #so 
    y += B
    
    return ylist, y


def chopperfunction(x_values, r_beam, w_slit):
    """x_values: if you want n amplitudes, length x_values is ceil(1/2n+1) and numbers range from -1/2w_slit to R"""
    
    edge = x_values >= r_beam - w_slit
    center = x_values < r_beam - w_slit
    
    y_values = np.zeros(len(x_values))
    
    def integral(r_beam, x):
        return x*np.sqrt(r_beam**2-x**2) + r_beam**2*np.arcsin(x/r_beam)
    
    y_values[edge] = integral(r_beam, r_beam) - integral(r_beam, x_values[edge])
    y_values[center] = (integral(r_beam, x_values[center]+w_slit) 
                        - integral(r_beam, x_values[center]))
    
    return y_values

def get_chopper_amplitudes(n=15, chopper_freq=250, r_beam=0.23, w_slit=0.85): #r and w in mm
    """ Returns only odd length array. if n is even, returns array of length n-1"""
    # xt_conversion = 24 * 250 / (w_slit * chopper_freq) #in microseconds, old conversion that is not correct
    xt_conversion = 1/(2*np.pi*55.992*chopper_freq) * 1e3 #1e3 to convert to milliseconds

    
    x_values = np.linspace(-w_slit/2, r_beam, num=int(np.ceil(n/2+1)))
    
    amplitudes = chopperfunction(x_values, r_beam, w_slit)
    amplitudes = np.concatenate((np.flip(amplitudes)[:-1],amplitudes))
    amplitudes = amplitudes[1:-1]
    amplitudes /= np.pi * r_beam**2 #to normalize
    amplitudes[np.isnan(amplitudes)] = 1 #very ugly, but it works..

    
    t_values = (x_values+w_slit/2) * xt_conversion
    t_values = np.concatenate(((-1*np.flip(t_values)[:-1]),t_values))
    t_values = t_values[1:-1]
    
    return t_values, amplitudes

main()