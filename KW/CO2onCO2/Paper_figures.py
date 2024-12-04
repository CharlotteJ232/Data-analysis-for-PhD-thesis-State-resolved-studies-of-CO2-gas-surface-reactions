import numpy as np
from matplotlib import pyplot as plt
import matplotlib.patches as patches
from matplotlib.ticker import (AutoMinorLocator,MultipleLocator)
import os


folderstart = "C:/Users/jansenc3/surfdrive/DATA/"
folderstart = "C:/Users/Werk/Surfdrive/DATA/"
# folder = '2022/01 Jan/220114/KW/' #old
folder = '2023/05 May/230508/Pressure/'
savefolder = 'CO2_on_CO2/Figures/2023/'
# kwname = 'KW03' #old
kwname = 'onresonance_R4'
lasername = 'onresonance_R4_laser'
# start = 75
# stop = 450
start = 1
stop = 5200
save_figures=False


def main():
    if save_figures:
        if not os.path.exists(folderstart+savefolder):
            os.makedirs(folderstart+savefolder)

    time, data = np.loadtxt(folderstart+folder+kwname+'.txt', skiprows=0, usecols=(0, 1), unpack=True)
    time -= time[0]
    time_laser, modulation = np.loadtxt(folderstart+folder+lasername+'.txt', skiprows=1, usecols=(0,1),unpack=True)
    time_laser -= time_laser[0]

    # kw_plot(time, data)
    # kw_plot(time, data, axis=[30,36,None,None])
    # kw_plot(time, data, axis=[474,482,None,None])
    
    #repeat normalization for both opening and closing flag 1
    norm1 = kw_normalize(time, data, low1=31,low2=33,high1=34,high2=36)
    norm2 = kw_normalize(time, data, low1=479,low2=481,high1=476,high2=478)
    norm = (norm1+norm2)/2
    #normalization placeholder until new measurement
    norm = 1.2E-8
    print(norm)

    #KW plot
    kw_plot(time, data/norm, start=start, stop=stop, save=save_figures, axis=[start,stop,0,1])

    #FFT plot
    fft_plot(xmin=4.9, xmax=5.1, save=save_figures, norm=norm)

    #mask plot
    mask_plot(time,time_laser, modulation, data, axis=[300, 301, None,None], save=save_figures)

    #difference plot
    difference_plot(norm=norm, axis=[-0.3,0.3, -0.0008,0.0008], save=save_figures)







#KW plot
def kw_plot(time, data, axis=[None,None,None,None], start=None, stop=None, save=False):
    fig,ax = plt.subplots()
    ax.plot(time, data)
    ax.axis(axis)

    #ticks and labels
    if axis[3]:
        plt.yticks(np.linspace(axis[2],axis[3],5))
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('CO2 partial pressure (normalized)')
    ax.tick_params(top=True, direction='in')  
    ax.tick_params(right=True, direction='in')
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(which='minor', top=True, direction='in')  
    ax.tick_params(which='minor', right=True, direction='in')


    if start and stop:
        rect = patches.Rectangle((start, 0.4), (stop-start), 0.7, linewidth=1, edgecolor='gray', facecolor='none')
        ax.add_patch(rect)

    if save:
        plt.savefig(folderstart+savefolder+'KW_plot.png', dpi=500)
    plt.show()
    plt.close()

#FFT plot
def fft_plot(xmin=None,xmax=None, save=False, norm=1, freq=3.0, freq_tolerance=1E-4):
    fftfreq_laser, fft_laser = np.loadtxt(folderstart+folder+kwname+'/fft_laser_'+str(start)+'_'+str(stop)+'.txt', unpack=True, skiprows=1)
    fftfreq, fft = np.loadtxt(folderstart+folder+kwname+'/fft_'+str(start)+'_'+str(stop)+'.txt', unpack=True, skiprows=1)
    #Check if fftfreq_laser and fftfreq are the same
    if np.any(fftfreq != fftfreq_laser):
        print ("warning! fftfreq and fftfreq_laser are not the same")
        diff = fftfreq != fftfreq_laser
        print (fftfreq[diff])
    #normalize fft with respect to KW measurement
    fft /= norm
    #normalize fft_laser with respect to modulation measurement
    fft_laser /= 5
    #normalize fft to reflect the amplitude of the sine 
    fft /= len(fftfreq)
    fft_laser /= len(fftfreq_laser)
    #now correct for the fact that the expected signal is a square wave, which amplitude is 
    #smaller than the amplitude of the sine wave at the same frequency
    fft *= np.pi/4
    fft_laser *= np.pi/4 

    #integration over the entire peak to get the actual value?
    fmin = 4.933 #chosen such that integration over fft_freq peak is exactly 5, which is the real modulation amplitude
    fmax = 5.067
    fmin = 2.933
    fmax = 3.067
    window = (fftfreq_laser > fmin) * (fftfreq_laser < fmax)
    print("laser peak integral stuff \n", np.sum(fft_laser[window]))

    # #integrate fft * fft_laser to integrate over the part of fft with the same shape as fft_freq
    # #something is not right here
    # print(np.sum(fft[window]*fft_laser[window])/np.sum(fft_laser[window]))

    # #maybe this is the right way to normalize
    # center = (fftfreq_laser > freq-0.001) * (fftfreq_laser < freq+0.001)
    # print(fft_laser[center])
    # print(np.sum(fft[window]*fft_laser[window])/np.sum(fft_laser[center]))

    # #or maybe this.. 
    center = (fftfreq_laser > freq-freq_tolerance) * (fftfreq_laser < freq+freq_tolerance)
    # print(fft[center])
    # print(fft[center]/fft_laser[center]*np.sum(fft_laser[window]))
    # print(np.sum(fft[center]*fft_laser[window])/np.sum(fft_laser[center]))

    
    #ACTUAL PLOT START
    fig, (ax1,ax2) = plt.subplots(2,1,sharex='all')
    fig.subplots_adjust(hspace=0)

    #plotted here is the absolute value of the FFT results
    ax1.plot(fftfreq, fft, label='Sticking probability FFT') 
    # ax1.fill_between(fftfreq[window],0,fft[window]*fft_laser[window]/np.sum(fft_laser[center]), color='lightgray')
    ax1.fill_between(fftfreq[window],0,fft[center]*fft_laser[window]/np.sum(fft_laser[center]), color='lightgray')
    # ax1.fill_between(fftfreq[window],0,fft[window]*fft_laser[window]/np.sum(fft_laser[window]), color='gray')
    ax2.plot(fftfreq_laser, fft_laser, label='Laser modulation FFT')
    ax2.fill_between(fftfreq[window],0,fft_laser[window], color='lightgray')

    ax1.text(0.025, 0.82, 'a)', transform=ax1.transAxes)
    ax2.text(0.025, 0.82, 'b)', transform=ax2.transAxes)

    ax2.set_xlabel('Frequency (Hz)')

    #axis limits
    ax1.set_ylim(0,0.0006)
    ax2.set_ylim(0, 0.6)
    ax2.set_xlim(xmin, xmax)
    ax2.tick_params(right=True, labelright=True, labelleft=False, direction='in')

    for ax in [ax1, ax2]:
        ax.plot([5,5],[0,100000],'--',c='gray')
        ax.tick_params(top=True, right=True, direction='in')  
        # ax.tick_params(right=False, left=False, labelleft=False, direction='in')
        ax.xaxis.set_minor_locator(AutoMinorLocator())
        ax.yaxis.set_minor_locator(AutoMinorLocator())
        ax.xaxis.set_major_locator(MultipleLocator((xmax-xmin)/4))
        ax.tick_params(which='minor', top=True, right=True, direction='in')  
        ax.legend()

    if save:
        plt.savefig(folderstart+savefolder+'fft_plot.png', dpi=500)
    plt.show()
    plt.close()

#mask plot
def mask_plot(time, time_laser, modulation, data, axis=[None,None,None,None], save=False):
    print('mask plot')
    fig, (ax1,ax2) = plt.subplots(2,1, sharex='all')
    fig.subplots_adjust(hspace=0)
    modulation_norm = np.max(modulation)
    ax1.plot(time, data, label='CO$_2$ partial pressure')
    ax2.plot(time_laser, modulation/modulation_norm, label='Laser modulation')

    ax1.text(0.025, 0.82, 'a)', transform=ax1.transAxes)
    ax2.text(0.025, 0.82, 'b)', transform=ax2.transAxes)

    #ticks and labels
    ax2.set_xlabel('Time (s)')
    # ax1.set_ylabel('CO2 partial pressure')
    # ax2.set_ylabel('Data mask')
    # ax3.set_ylabel('Shifted data mask')
    for ax in [ax1, ax2]:
        ax.tick_params(top=True, direction='in')  
        ax.tick_params(left=False, labelleft=False, direction='in')
        ax.xaxis.set_minor_locator(AutoMinorLocator())
        ax.tick_params(which='minor', top=True, direction='in')  
        ax.legend(loc='lower right')

    #axis limits
    ax1.set_ylim(5.65E-9,5.75E-9)
    ax2.set_ylim(-0.5, 1.5)
    ax2.set_xlim(axis[0],axis[1])

    if save:
        plt.savefig(folderstart+savefolder+'mask_plot.png', dpi=500)
    plt.show()
    plt.close()

#mask plot 2
def mask_plot2(time, modulation, data, axis=[None,None,None,None], save=False):
    print('mask plot 2')
    fig, (ax1,ax2,ax3) = plt.subplots(3,1, sharex='all')
    fig.subplots_adjust(hspace=0)
    modulation_norm = np.max(modulation)
    ax1.plot(time, data, label='Sticking probability data')
    ax2.plot(time, modulation/modulation_norm, label='Data mask')
    ax3.plot(time+0.05, modulation/modulation_norm, label='Data mask (shifted by 0.05 s)')

    #ticks and labels
    ax3.set_xlabel('Time (s)')
    # ax1.set_ylabel('CO2 partial pressure')
    # ax2.set_ylabel('Data mask')
    # ax3.set_ylabel('Shifted data mask')
    for ax in [ax1, ax2, ax3]:
        ax.tick_params(top=True, direction='in')  
        ax.tick_params(left=False, labelleft=False, direction='in')
        ax.xaxis.set_minor_locator(AutoMinorLocator())
        ax.tick_params(which='minor', top=True, direction='in')  
        ax.legend(loc='lower right')

    #axis limits
    ax1.set_ylim(6.2,6.5)
    ax2.set_ylim(-0.2, 1.2)
    ax3.set_ylim(-0.2, 1.2)
    ax3.set_xlim(axis[0],axis[1])

    if save:
        plt.savefig(folderstart+savefolder+'mask_plot2.png', dpi=500)
    plt.show()
    plt.close()

#difference plot
def difference_plot(norm=1, axis=[None,None,None,None], save=False, freq=5):
    shift, difference = np.loadtxt(folderstart+folder+kwname+'/timeshift_'+str(start)+'_'+str(stop)+'.txt', unpack=True, skiprows=1)
    fig, ax = plt.subplots()
    ax.plot(shift, difference/norm)
    ax.axis(axis)

    ax.plot([-5,5],[0,0],'--',c='gray')

    #ticks and labels
    if axis[2] and axis[3]:
        plt.yticks(np.linspace(axis[2],axis[3],5))
    ax.set_xlabel('Time shift of the data mask (s)')
    ax.set_ylabel('Measured sticking probability difference')
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(top=True, direction='in')  

    ax.tick_params(which='minor', top=True, direction='in')  
    # ax.tick_params(which='minor', right=True, direction='in')

    #axis with absolute differences corrected for J=2 population etc
    population = 0.2
    excitation = 1
    factor = population*excitation
    ax2 = ax.twinx()
    y2min = axis[2]/factor
    y2max = axis[3]/factor
    ax2.set_ylim(y2min, y2max)
    ax2.yaxis.set_major_locator(MultipleLocator((y2max-y2min)/4))
    ax2.yaxis.set_minor_locator(AutoMinorLocator())
    ax2.set_ylabel('Corrected sticking probability difference')
    ax2.tick_params(right=True, direction='in')
    ax2.tick_params(which='minor', right=True, direction='in')




    if save:
        plt.savefig(folderstart+savefolder+'difference_plot.png', dpi=500)
    plt.show()
    plt.close()


def kw_normalize(time, data, low1=None, low2=None, high1=None, high2=None):
    """
    low1, low2, time bounds for where the pressure is low (a small area on one side of the pressure jump)
    high1, high2 time bounds for where the pressure is high (close to the low bounds)
    """
    masklow = (time>low1) * (time<low2)
    maskhigh = (time>high1) * (time<high2)
    return np.average(data[maskhigh]) - np.average(data[masklow])




main()