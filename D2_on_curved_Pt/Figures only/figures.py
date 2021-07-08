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




def main():
    slope_and_111()
    s0_vs_pos()
    


def s0_vs_pos_img():
    fig, (ax_im, ax) = plt.subplots(2,1, figsize=(6,6))
    fig.subplots_adjust(hspace=0)
    position, position_centered, dens = np.loadtxt(savefolder+'step_density.txt', skiprows=1,unpack=True)
    
    colorcounter = 0
    for energy in list(datadic.keys()):
        color = cm.CMRmap(colorcounter/len(list(datadic.keys())))
        colorcounter += 1
        # color = cc.CET_R1[int(energy/np.max(list(datadic.keys()))*(len(cc.CET_R1)-1)*1)]
        # color = cm.gnuplot(energy/np.max(list(datadic.keys())))
        pos, s0, error = np.loadtxt(folderstart+year+datadic[energy]['kw_day']+file_kw, unpack=True, skiprows=1)
        ax.errorbar(dens, s0, error, label=str(energy), color=color,fmt='o',elinewidth=5)
     
    def pos_to_dens(x):
        return np.interp(x, position_centered,dens)
        
    def dens_to_pos(x):
        return np.interp(x, dens,position_centered)

    #axis manipulation graph
    ax.tick_params(top=False)
    ax.set_xlim(-0.5,0.5) #cannot be higher unless I calculate a broader dens_to_pos for interpolation
    ax.set_xlabel('Step density (/nm)')
    ax.set_ylabel('s0')
    ax.locator_params(axis='y',nbins=5)
    ax.locator_params(axis='x', nbins=6)
    ax.set_ylim(0,0.4)
    densticks = ax.get_xticks()
    ax.set_xticklabels(np.round(np.absolute(densticks),1)) #to get absolute step densities on both sides
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
   
    #image

    
   #secondary axis with crystal position 
    xmin, xmax = ax.get_xlim()
    ax2=ax_im.twiny()
    ax_im.set_xlim(xmin, xmax)
    xmin,xmax=ax.get_xlim()
    x2min=dens_to_pos(xmin)
    x2max=dens_to_pos(xmax)
    ax2.set_xlim(x2min,x2max)
    ax2.set_xlabel('Position (mm)')
    ax2.xaxis.set_major_locator(MaxNLocator(nbins=6))
    ax2.xaxis.set_minor_locator(AutoMinorLocator())
    #other things
    ax.legend(bbox_to_anchor=(1.02, 0,1.02,1), borderaxespad=0,loc=2)
    
    
    img = mpimg.imread(image_path)
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    ax_im.axis([xmin, xmax, ymin,ymax])
    ax2.imshow(img, extent=[x2min, x2max, ymin, ymax])
    ax_im.axis([xmin, xmax, ymin,ymax])
    ax_im.tick_params(bottom=False, top=False, left=False, right=False, labelleft=False, labelbottom=False)

    
    plt.show()
    plt.close()    
    
def s0_vs_pos():
    fig,  ax = plt.subplots()

    position, position_centered, dens = np.loadtxt(savefolder+'step_density.txt', skiprows=1,unpack=True)
    

    
    colorcounter = 0
    for energy in list(datadic.keys()):
        color = cm.CMRmap(colorcounter/len(list(datadic.keys())))
        colorcounter += 1
        # color = cc.CET_R1[int(energy/np.max(list(datadic.keys()))*(len(cc.CET_R1)-1)*1)]
        # color = cm.gnuplot(energy/np.max(list(datadic.keys())))
        pos, s0, error = np.loadtxt(folderstart+year+datadic[energy]['kw_day']+file_kw, unpack=True, skiprows=1)
        ax.errorbar(dens, s0, error, label=str(energy), color=color,fmt='o',elinewidth=5)
     
    def pos_to_dens(x):
        return np.interp(x, position_centered,dens)
        
    def dens_to_pos(x):
        return np.interp(x, dens,position_centered)

    #axis manipulation
    ax.tick_params(top=False)
    ax.set_xlim(-0.5,0.5) #cannot be higher unless I calculate a broader dens_to_pos for interpolation
    ax.set_xlabel('Step density (/nm)')
    ax.set_ylabel('s0')
    ax.locator_params(axis='y',nbins=5)
    ax.locator_params(axis='x', nbins=6)
    ax.set_ylim(0,0.4)
    densticks = ax.get_xticks()
    ax.set_xticklabels(np.round(np.absolute(densticks),1)) #to get absolute step densities on both sides
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    #secondary axis with crystal position
    ax2=ax.twiny()
    xmin,xmax=ax.get_xlim()
    x2min=dens_to_pos(xmin)
    x2max=dens_to_pos(xmax)
    ax2.set_xlim(x2min,x2max)
    ax2.set_xlabel('Position (mm)')
    ax2.xaxis.set_major_locator(MaxNLocator(nbins=6))
    ax2.xaxis.set_minor_locator(AutoMinorLocator())
    #other things
    ax.legend(bbox_to_anchor=(1.02, 0,1.02,1), borderaxespad=0,loc=2)

    
    plt.show()
    plt.close()
    
def slope_and_111():
    energy, slopeB, slopeA, s111 = np.loadtxt(fitfolder+'slope_and_111.txt', skiprows=1, unpack=True)
    
    #plot slope and s111
    fig, (ax1, ax2, ax3) = plt.subplots(3,1, sharex='col',figsize=(6,6))
    fig.subplots_adjust(hspace=0)
    ax1.scatter(energy, slopeA, color='blue', label='A type')
    ax1.scatter(energy,slopeB, color='red', label='B type')
    ax2.scatter(energy, s111, color='black', label='s0 (111)')
    
    ax1.legend(loc='upper right')
    ax2.legend(loc='lower right')
    
    ax1.axis([0,8,0,0.2])
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


main()