# -*- coding: utf-8 -*-
"""
https://matplotlib.org/3.2.1/tutorials/introductory/customizing.html
https://stackoverflow.com/questions/26106552/matplotlib-style-library-not-updating-when-mplstyle-files-added-deleted
https://matplotlib.org/3.2.1/tutorials/introductory/customizing.html#customizing-with-matplotlibrc-files
@author: Charlotte
"""
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.ticker import *
from matplotlib import cm
import matplotlib.image as mpimg
import colorcet as cc


plt.style.reload_library()
plt.style.use('voorbeeld')

folderstart = 'P:/Surfdrive/DATA/'
folderstart = 'D:/Surfdrive/DATA/'
year = '2020/'
file_kw = '/KW/Images/stickingprob_vs_position.txt'
file_TOF = '/TOF/Images/4.0/energy.txt'
savefolder = folderstart+'D2 sticking probability as a function of step density and kinetic energy/Images/'
fitfolder = savefolder
image_path = "D:/Surfdrive/DATA/D2 sticking probability as a function of step density and kinetic energy/Images/curved_crystal_placeholder.png"

#File names and locations for s0 data
datadic = {0.7:{'kw_day':'03 Mar/200310', 'TOF_day':'03 Mar/200310','c':'red'},
           0.8:{'kw_day':'05 May/200512', 'TOF_day':'01 Jan/200128','c':'orange'},
           1.2:{'kw_day':'03 Mar/200305', 'TOF_day':'02 Feb/200220','c':'yellow'},
           2.0:{'kw_day':'06 Jun/200602', 'TOF_day':'05 May/200511','c':'green'},
           2.9:{'kw_day':'03 Mar/200303', 'TOF_day':'02 Feb/200218','c':'blue'},
           4.0:{'kw_day':'03 Mar/200311', 'TOF_day':'03 Mar/200306','c':'indigo'},
           4.77:{'kw_day':'05 May/200519', 'TOF_day':'05 May/200519','c':'purple'},
           6.1:{'kw_day':'05 May/200520', 'TOF_day':'02 Feb/200204','c':'violet'}}

#fitted slopes and s111 data


ext='.pdf'

def main():
    plot_s0_vs_E(datalist)
    


def plot_s0_vs_E(datalist, use_stepdens=True, save=False): 
    """
    This function:
        - fits an exponential + line through the s0 vs E data
        - plots the data + fit
        - adds the fitted s0t to the data objects, as an array with elements for each position/step density
    """
    if use_stepdens:
        dic = {k:[] for k in datalist[0].probed_step_density}
    else:
        dic = {k:[] for k in datalist[0].pos}
    E = []
    for dataset in datalist:
        dataset.s0t = [] #make empty list for s0t (vs position), will be filled after fitting
        E.append(dataset.E_avg)
        for i in range(len(dataset.data)):
            if use_stepdens:
                dic[dataset.probed_step_density[i]].append(dataset.data[i])
            else:
                dic[dataset.pos[i]].append(dataset.data[i])  
    
    #FIT AND PLOT
    
    #for colormap
    stepdensities = np.array(list(dic.keys()))
    max_density = np.max(np.absolute(stepdensities))
    def to_color(step_density):
        return cm.RdBu(step_density/max_density*0.5+0.5)
    
    #plot
    fig, axes = plt.subplots(2,1,sharex='col',figsize=(10,10))
    ax1 = axes[0]
    ax2 = axes[1]
    fig.subplots_adjust(hspace=0)
    
    #FIT
    E = np.array(E)
    Efit = np.arange(np.min(E),np.max(E),0.01)
    def expline(E, C, tau, a, b):
        return C*np.exp(-tau*E)+a*E+b
        
    guess = [0.5, 100, 1, 0.15] #values for E in eV, taken from groot 2011 supporting information
    guess = [0.5, 1, 0.01, 0.15] #values for E in kJ/mol
    lower = [-np.inf, -np.inf, -np.inf, 0]
    upper = [np.inf, np.inf, np.inf, np.inf]
    
    for k in dic.keys():
        #FIT
        popt, pcov = curve_fit(expline, E, dic[k],p0=guess, bounds=(lower,upper))
        #ADD S0T TO DATA OBJECTS
        for dataset in datalist:
            dataset.s0t.append(popt[2]*dataset.E_avg)
        
        #PLOT ORIGINAL AND FIT
        if k<=0:
            # ax2.plot(E,popt[2]*E+popt[3],color=to_color(k)) #plot just the line
            ax2.plot(E, dic[k],'o', color=to_color(k)) 
            ax2.plot(Efit, expline(Efit,*popt),'-',label=str(k), color=to_color(k))
        else:
            ax1.plot(E, dic[k],'o',label=str(k), color=to_color(k))
            ax1.plot(Efit, expline(Efit,*popt),'-',label=str(k), color=to_color(k))
     
    for dataset in datalist:
        dataset.s0t=np.array(dataset.s0t)        
     
    #axis limits
    axlims = [0,None,0,0.45]
    ax1.axis(axlims)
    ax2.axis(axlims)    
    
    #make colorbar
    sm = plt.cm.ScalarMappable(cmap=cm.RdBu, norm=plt.Normalize(vmin=-max_density, vmax=max_density))
    cbar = fig.colorbar(sm, ax=axes.ravel().tolist())
    labels = cbar.ax.get_yticks()
    cbar.ax.set_yticklabels(np.round(np.absolute(labels),1))
    cbar.ax.set_ylabel('step density (/nm)')   

    #labels and title
    ax2.set_xlabel('Energy (kJ/mol)')
    ax1.set_ylabel('s0')
    ax2.set_ylabel('s0')
    ax1.set_title('B type (red) and A type (blue)')
    if save:
        plt.savefig(savefolder+'s0_vs_E.png',dpi=500)
    plt.show()        

    
    
    
def slope_and_111(fitranges=1.2, linestyle='-o', errobars=False, save=False):
    
    #plot slope and s111
    fig, (ax1, ax2, ax3) = plt.subplots(3,1, sharex='col',figsize=(6,6))
    fig.subplots_adjust(hspace=0)
    
    if errobars:
        A_list = []
        B_list = []
        s111_list = []
        for fitrange in fitranges:
            energy, slopeB, slopeA, s111 = np.loadtxt(fitfolder+'slope_and_111_range_'+str(fitrange)+'.txt', skiprows=1, unpack=True)
            A_list.append(slopeA)
            B_list.append(slopeB)
            s111_list.append(s111)
        slopeA = np.mean(A_list, axis=0)
        ErrorA = np.std(A_list, axis=0)
        slopeB = np.mean(B_list, axis=0)
        ErrorB = np.std(B_list, axis=0)
        s111 = np.mean(s111_list, axis=0)
        Error111 = np.std(s111_list, axis=0)        
        ax1.errorbar(energy, slopeA, ErrorA, color='blue', fmt='o', label='A type')
        ax1.errorbar(energy, slopeB, ErrorB, color='red', fmt='o', label='B type')
        ax2.errorbar(energy, s111, Error111, color='black', fmt='o', label='s0 (111)')            
    
    else:
        for fitrange in fitranges:
            energy, slopeB, slopeA, s111 = np.loadtxt(fitfolder+'slope_and_111_range_'+str(fitrange)+'.txt', skiprows=1, unpack=True)
        
            ax1.plot(energy, slopeA, linestyle, color='blue', label='A type')
            ax1.plot(energy,slopeB, linestyle, color='red', label='B type')
            ax2.plot(energy, s111, linestyle, color='black', label='s0 (111)')
    
    ax1.legend(loc='upper right')
    ax2.legend(loc='lower right')
    
    ax1.axis([0,8,0,0.3])
    ax2.axis([None,None,0,0.2])
    
    ax1.tick_params(top=False)
    ax2.tick_params(labelleft=False, labelright=True)

    ax1.set_ylabel('Slope')
    ax2.set_ylabel('s0')
    ax2.yaxis.set_label_position('right')
    

    
    #Plot energy 
    for energy in list(datadic.keys()):
        E, y = np.loadtxt(folderstart+year+datadic[energy]['TOF_day']+file_TOF, skiprows=1, unpack=True)
        ax3.plot(E,y, color=cm.gray(energy/10))
        ax3.scatter(energy, np.interp(energy, E,y),color=cm.gray(energy/10))
    ax3.plot(-1,-1,color='grey',label='Energy distribution')
    ax3.scatter(-1,-1,color='grey',label='Average energy')
    ax3.set_ylim(0,0.01)
    ax3.set_yticklabels([])
    ax3.set_ylabel('Probability density')
    ax3.legend(loc='upper right')
    ax3.set_xlabel('Energy (kJ/mol)')
    

    
    def kjtoev(x):
        return x/96.485*1000
    
    def evtokj(x):
        return x/1000*96.485
    
    secax1 = ax1.secondary_xaxis('top', functions=(kjtoev, evtokj))
    secax1.set_xlabel('Energy (meV)')
    
    fig.show()
    
    if save:
        fig.savefig(savefolder+'fig_slopes_111'+ext)


main()