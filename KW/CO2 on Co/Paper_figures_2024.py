import numpy as np
from matplotlib import pyplot as plt
import matplotlib.patches as patches
from matplotlib.ticker import (AutoMinorLocator,MultipleLocator)
import os
import matplotlib as mpl


folderstart = "C:/Users/jansenc3/surfdrive/"
folderstart = "C:/Users/Werk/Surfdrive/"
folder = 'DATA/2024/02 Feb/240222/kw/'
savefolder = 'Proefschrift/CO2_on_CO2_physisorption/Figures/'
savefolder = 'Proefschrift/Cobalt/Figures/'
# savefolder = folder + 'Figures/'
kwname = 'KW01'
lasername = 'KW01_laser'
filename_normalize = 'KW01'
folder_normalize = folder
start = 160
stop = 6300
save_figures=True

plotlines = [(136.5, 'open flag 1'), (153, 'open flag 2'), (6456.5, 'close flag 2'), (6475.5, 'close flag 1')]
vertical_duration = 1.5 #s
average_over = 0.5 #s




def main():
    if save_figures:
        if not os.path.exists(folderstart+savefolder):
            os.makedirs(folderstart+savefolder)

    time, data = np.loadtxt(folderstart+folder+kwname+'.txt', skiprows=3, usecols=(0, 1), unpack=True)
    time -= time[0]
    time_laser, modulation = np.loadtxt(folderstart+folder+lasername+'.txt', skiprows=1, usecols=(0,1),unpack=True)
    time_laser -= time_laser[0]
    time_normalize, data_normalize = np.loadtxt(folderstart+folder_normalize+filename_normalize+'.txt', skiprows=3, usecols=(0,1), unpack=True)
    time_normalize -= time_normalize[0]


    kw_plot(time_normalize, data_normalize)
    # kw_plot(time_normalize, data_normalize, axis=[20,25,None,None])
    # kw_plot(time_normalize, data_normalize, axis=[120,125,None,None])
    # kw_plot(time, data, axis=[474,482,None,None])
    
    #repeat normalization for both opening and closing flag 1

    norms = []
    mask = np.full(len(time_normalize), False)
    masks = []
    for line in plotlines:
        tim = line[0]
        norm_temp, mask_temp = kw_normalize(time_normalize, data_normalize, low1=tim-average_over, low2=tim, high1=tim+vertical_duration, high2=tim+vertical_duration+average_over)
        norms.append(norm_temp)
        mask = mask | mask_temp
        masks.append(mask_temp)

    print('flag1, flag2, flag2, flag1 ',norms)
    norm = (norms[0]+norms[3])/2 #calculate normalisation for plot
    S0 = [norms[1]/norms[0], norms[2]/norms[3]]
    print('S (begin and end, not real though) ', S0)

    # #KW plot
    # set_zero = np.min(data_normalize[0:int(len(data_normalize)/3)])
    # kw_plot(time_normalize, (data_normalize-set_zero)/norm+0.1, start=0, stop=150, save=save_figures, axis=[0,150,0,1.6], plotlines=plotlines, plotmasks=masks)
    # # kw_plot(time, data/norm, start=start, stop=stop, save=save_figures, axis=[start,stop,0,1])

    #FFT plot
    s_max = fft_plot(xmin=2.986, xmax=3.014, save=save_figures, norm=norm)

    rot_pop = 0.2 #rotational state population
    emission = 0.7 #population accounted for spontaneous emission
    s_max_corrected = s_max / rot_pop / emission
    print('actual s_max', s_max_corrected)



    # #mask plot
    # mask_plot(time,time_laser, modulation, data, axis=[300, 301, None,None], save=save_figures)

    # # difference plot
    # difference_plot(norm=norm, axis=[-0.3,0.3, -0.0008,0.0008], save=save_figures)







#KW plot
def kw_plot(time, data, axis=[None,None,None,None], start=None, stop=None, save=False, plotlines=[], plotmasks=None):
    fig,ax = plt.subplots()
    
    if axis[3]:
        ytext = axis[3]*1.02
        for lineinfo in plotlines:
            plot_line(ax, *lineinfo, ytext, fontsize=10)
    
    ax.plot(time, data)
    if plotmasks is not None:
        for mask in plotmasks:
            ax.plot(time[mask], data[mask], color='aqua')
    ax.axis(axis)

    #ticks and labels
    if axis[3]:
        plt.yticks(np.linspace(axis[2],axis[3],5))
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('CO$_2$ partial pressure (normalized)')
    ax.tick_params(top=True, direction='in')  
    ax.tick_params(right=True, direction='in')
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(which='minor', top=True, direction='in')  
    ax.tick_params(which='minor', right=True, direction='in')
    figure_look(ax)



    if start and stop:
        rect = patches.Rectangle((start, 0.4), (stop-start), 0.7, linewidth=1, edgecolor='gray', facecolor='none')
        ax.add_patch(rect)


    
    if save:
        plt.savefig(folderstart+savefolder+'KW_plot.png', bbox_inches='tight', dpi=500)
        plt.savefig(folderstart+savefolder+'KW_plot.eps', bbox_inches='tight', dpi=500)

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
    fmin = 4.933 #chosen such that integration over fft_freq peak is exactly 1? or 5?, which is the real modulation amplitude
    fmax = 5.067
    fmin = 2.999
    fmax = 3.001
    window = (fftfreq_laser > fmin) * (fftfreq_laser < fmax)
    print("laser peak integral stuff \n", 'laser peak integral (should be 1)', np.sum(fft_laser[window]))

    # #integrate fft * fft_laser to integrate over the part of fft with the same shape as fft_freq
    # #something is not right here
    # print(np.sum(fft[window]*fft_laser[window])/np.sum(fft_laser[window]))

    # #maybe this is the right way to normalize
    # center = (fftfreq_laser > freq-0.001) * (fftfreq_laser < freq+0.001)
    # print(fft_laser[center])
    # print(np.sum(fft[window]*fft_laser[window])/np.sum(fft_laser[center]))

    # #or maybe this.. 
    # center = (fftfreq_laser > freq-freq_tolerance) * (fftfreq_laser < freq+freq_tolerance)
    # print('fft[center]',fft[center])
    # print('fft_center normalized', fft[center]/fft_laser[center]*np.sum(fft_laser[window]))
    # print(np.sum(fft[center]*fft_laser[window])/np.sum(fft_laser[center]))

    #or simply this
    center = (fftfreq_laser > freq-freq_tolerance) * (fftfreq_laser < freq+freq_tolerance)
    s_max = fft[center]/fft_laser[center]
    print('this is the maximum change in sticking probability (fft center/fft laser center): ', s_max)

    
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

    ax1.text(0.025, 0.83, 'a', transform=ax1.transAxes)
    ax2.text(0.025, 0.83, 'b', transform=ax2.transAxes)

    ax2.set_xlabel('Frequency (Hz)')

    #axis limits
    ax1.set_ylim(0,2E-4)
    ax2.set_ylim(0, 0.6)
    ax2.set_xlim(xmin, xmax)
    ax2.tick_params(right=True, labelright=True, labelleft=False, direction='in')

    for ax in [ax1, ax2]:
        figure_look(ax)
        ax.plot([5,5],[0,100000],'--',c='gray')
        ax.tick_params(top=True, right=True, direction='in')  
        # ax.tick_params(right=False, left=False, labelleft=False, direction='in')
        ax.xaxis.set_minor_locator(AutoMinorLocator())
        ax.yaxis.set_minor_locator(AutoMinorLocator())
        ax.xaxis.set_major_locator(MultipleLocator((xmax-xmin)/4))
        ax.tick_params(which='minor', top=True, right=True, direction='in')  
        ax.legend(fontsize=10, loc='upper right')

    if save:
        plt.savefig(folderstart+savefolder+'fft_plot.png', dpi=500)
        plt.savefig(folderstart+savefolder+'fft_plot.eps', dpi=500)
        plt.savefig(folderstart+savefolder+'fft_plot.pdf', dpi=500)
    plt.show()
    plt.close()

    return s_max

#mask plot
def mask_plot(time, time_laser, modulation, data, axis=[None,None,None,None], save=False):
    print('mask plot')
    fig, (ax1,ax2) = plt.subplots(2,1, sharex='all')
    fig.subplots_adjust(hspace=0)
    modulation_norm = np.max(modulation)

    ax1.plot(time-axis[0], data, label='CO$_2$ partial pressure')
    ax2.plot(time_laser-axis[0], modulation/modulation_norm, label='Laser modulation')

    # ax1.text(0.025, 0.82, 'a)', transform=ax1.transAxes)
    # ax2.text(0.025, 0.82, 'b)', transform=ax2.transAxes)

    #ticks and labels
    ax2.set_xlabel('Time (s)')
    # ax1.set_ylabel('CO2 partial pressure')
    # ax2.set_ylabel('Data mask')
    # ax3.set_ylabel('Shifted data mask')
    for ax in [ax1, ax2]:
        figure_look(ax)
        ax.tick_params(top=True, direction='in')  
        ax.tick_params(left=False, labelleft=False, direction='in')
        ax.xaxis.set_minor_locator(AutoMinorLocator())
        ax.tick_params(which='minor', top=True, direction='in')  
        ax.legend(loc='lower right', fontsize=10)

    #axis limits
    ax1.set_ylim(5.65E-9,5.75E-9)
    ax2.set_ylim(-0.5, 1.5)
    ax2.set_xlim(0,axis[1]-axis[0])

    if save:
        plt.savefig(folderstart+savefolder+'mask_plot.png', dpi=500)
        plt.savefig(folderstart+savefolder+'mask_plot.eps', dpi=500)
    plt.show()
    plt.close()



def plot_line(ax, wavenumber, text, ytext, fontsize=14):
    ax.plot([wavenumber,wavenumber],[-100,100], '--', lw=1, color='gray')
    ax.text(wavenumber, ytext, text, size=fontsize, rotation=45)


def kw_normalize(time, data, low1=None, low2=None, high1=None, high2=None):
    """
    low1, low2, time bounds for where the pressure is low (a small area on one side of the pressure jump)
    high1, high2 time bounds for where the pressure is high (close to the low bounds)
    """
    masklow = (time>low1) * (time<low2)
    maskhigh = (time>high1) * (time<high2)
    if low1<high1:
        masktotal = (time>((low1+low2)/2))*(time<((high1+high2)/2))
    else: #this is not necessary anymore as low and high are meaningless when I use np.absolute anyway
        masktotal = (time<((low1+low2)/2))*(time>((high1+high2)/2))
    return np.absolute(np.average(data[maskhigh]) - np.average(data[masklow])), np.array(masktotal)

def figure_look(ax, fontsize=12):
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(axis='x', top=True, direction='in')  #top=True, 
    ax.tick_params(axis='y', left=True, right=True, direction='in') #, labelleft=True, right=True, 
    ax.tick_params(which='minor', top=True, direction='in')  
    ax.tick_params(which='minor', left=True, right=True, direction='in')
    mpl.rcParams.update({'font.size': fontsize})




main()