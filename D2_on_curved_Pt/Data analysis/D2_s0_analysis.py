# -*- coding: utf-8 -*-
"""
Created on Fri Feb 28 15:05:09 2020

@author: jansenc3

NOTE: THIS DOES NOT HAVE THE FIX FOR FITRANGE PROBLEMS, COPY/PASTE THAT FROM 
THE COS VERSION OF THESCRIPT LATER
"""

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import cm
from scipy.optimize import least_squares, curve_fit
from matplotlib.offsetbox import TextArea, DrawingArea, OffsetImage, AnnotationBbox
from matplotlib.patches import Wedge, Circle, Arrow, Rectangle
import matplotlib.image as mpimg
import shutil
plt.style.reload_library()
#plt.style.use('voorbeeld')

unit_cell_width = 2*0.13874 #Pt atom diameter

folderstart = 'P:/Surfdrive/DATA/'
folderstart = 'C:/Users/Werk/Surfdrive/DATA/'
year = '2020/'
file_kw = '/KW/Images/stickingprob_vs_position.txt'
file_TOF = '/TOF/Images/4.0/corrected/energy.txt'
loadfolder = folderstart+'D2 sticking probability as a function of step density and kinetic energy/Images/'
savefolder = folderstart+'D2 sticking probability as a function of step density and kinetic energy/Images/Corrected/'
save_figs = False

datadic = {0.7:{'kw_day':'03 Mar/200310', 'TOF_day':'03 Mar/200310','c':'red'},
           0.8:{'kw_day':'05 May/200512', 'TOF_day':'01 Jan/200128','c':'orange'},
           1.2:{'kw_day':'03 Mar/200305', 'TOF_day':'02 Feb/200220','c':'gold'},
           2.0:{'kw_day':'06 Jun/200602', 'TOF_day':'05 May/200511','c':'yellow'},
           2.9:{'kw_day':'03 Mar/200303', 'TOF_day':'02 Feb/200218','c':'greenyellow'},
           4.0:{'kw_day':'03 Mar/200311', 'TOF_day':'03 Mar/200306','c':'green'},
           4.77:{'kw_day':'05 May/200519', 'TOF_day':'05 May/200519','c':'mediumseagreen'},
           6.1:{'kw_day':'05 May/200520', 'TOF_day':'02 Feb/200204','c':'teal'},
           7.8:{'kw_day':'08 Aug/200817', 'TOF_day':'08 Aug/200813_2','c':'blue'},
           9.4:{'kw_day':'08 Aug/200814', 'TOF_day':'08 Aug/200813','c':'darkblue'},
           10.7:{'kw_day':'07 Jul/200727', 'TOF_day':'07 Jul/200727','c':'indigo'},
           12.9:{'kw_day':'09 Sep/200904', 'TOF_day':'09 Sep/200904','c':'purple'},
           13.1:{'kw_day':'09 Sep/200910', 'TOF_day':'09 Sep/200910','c':'black'}}

datalist=[] #will be filled with data objects


fitrange = 1.0 #program fails for some values due to arrays not being the same size

def main():
    calc_attr()

    #read data and make data objects
    for dataset in datadic.values():
        datalist.append(data(folderstart, year, dataset['kw_day'],dataset['TOF_day'], 
                             file_kw, file_TOF, dataset['c'])) #list of data objects
    
    # # plot the energy distribution obtained by TOF measurements
    # plot_energy_dist(datalist)

    #test convolution
    test_convolution(datalist,save=save_figs)
    
    #plot the raw s0 data
    plot_s0_vs_pos(datalist,save=save_figs)

    #Fit straight lines through the data, now only mostly useful to fit the center
    center,slopeleft,sloperight,offset = fit_all_slopes(datalist,fitrange=fitrange,center_guess=19.91,save=False, plot=False)
    
    #make the dataset.probed_step_density array
    for dataset in datalist:
        dataset.convert_step_density(center, plot=False, save=False)
    
    
    #FIGURE 2
    #plot s0 vs E and fit the exp + line. necessary for data analysis. don't comment out 
    plot_s0_vs_E(datalist, use_stepdens=True, save=save_figs) 
    
    # #calculate cross section and calculate BA ratio, plot BA ratio
    # fig, ax = plt.subplots()
    # for dataset in datalist:
    #     dataset.calc_cross_section()
    #     dataset.save_crossec_s0(folderstart+'D2 sticking probability as a function of step density and kinetic energy/Crossec/')
    #     dataset.calc_BAratio()
    #     print(dataset.E_avg, dataset.BAratio)
    #     ax.scatter(dataset.E_avg,dataset.BAratio,color='black')
    # ax.plot(-1,-1,color='black',ls='none',marker='o',label='This study') #for legend only
    # ax.plot(0.9, 1.45,color='gray', marker='s',fillstyle='none',ls='none',label='van Lent et al.')
    # ax.plot(0.9, 1.54,color='gray', marker='v',fillstyle='none',ls='none',label='Auras et al.')
    # ax.plot([0,14],[1,1],ls='--',color='black')
    # ax.axis([0,14,0,2.5])
    # ax.tick_params(top=True, direction='in')  
    # ax.tick_params(right=True, direction='in')
    # plt.legend(loc='lower right')
    # # plt.title('BAratio')
    # plt.ylabel('B-type/A-type cross section ratio')
    # plt.xlabel('Energy (kJ/mol)')
    # if save_figs:
    #     plt.savefig(savefolder+'BAratio.png',dpi=500)
    #     plt.savefig(savefolder+'BAratio.eps')
    # plt.show()
    # plt.close()

    # #FIGURE 3
    # #plot cross section for each position individually
    # plot_cross_section_E(datalist, save=save_figs)

    # # #FIGURE 4
    # # # plot cross section as a function of step density
    # plot_cross_section_stepdensity(datalist, save=save_figs)
   
    #  # #SUPPLEMENTARY
    # plot_individual_cross_section(datalist, save=save_figs)
   
    # #SUPPLEMENTARY
    # #plot temperature dependence
    # plot_vs_temperature(save=save_figs)

    # # SUPPLEMENTARY
    # # hard cube comparison
    # plot_cross_section_E_hardcube(datalist,[6,12], save=save_figs)

    # #SUPPLEMENTARY
    # #compare to irenes cross section
    # compare_irene(datalist,save=save_figs)
    
    
"""
DATA CLASS 
"""    
class data:
    def __init__(self,folderstart,year,kw_day, TOF_day, file_kw, file_TOF, color):
        self.color=color
        self.pos, self.data, self.error= np.loadtxt(folderstart+year+kw_day+file_kw, 
                                                    unpack=True, usecols=(0,1,2), skiprows=1)
        if np.average(self.pos)>100:
            self.pos /= 10 #to convert to mm in some datasets
        self.energies, self.energy_dist = np.loadtxt(folderstart+year+TOF_day+file_TOF, 
                                                    unpack=True, usecols=(0,1), skiprows=1)
        self.E_avg = np.sum(self.energies*self.energy_dist)


    
    def plot_TOF(self):
        plt.plot(self.energies, self.energy_dist, c=self.color)
        plt.scatter(self.E_avg, np.interp(self.E_avg, self.energies, self.energy_dist),
                    c=self.color, label='avg = '+str(np.round(self.E_avg,decimals=2)))
    
    def plot_s0(self):
        plt.scatter(self.pos,self.data, color=self.color,label='E_avg = '+str(np.round(self.E_avg,decimals=2)))
        plt.errorbar(self.pos,self.data,yerr=self.error,capsize=5,color=self.color)#,fmt='none'
    
        

    def convert_step_density(self, center,save=False, plot=False):
        """
        Specific for beam of 124 um, assumes square beam profile and curved pt111 crystal
        Also calculates the terrace/step ratio (assumes single atom row for step)
        """
        step_height = 0.2265637926 #nm, for steps on pt111
        radius = 0.13874 #nm, for platinum atom
        crystal_curvature_radius = 15 #mm
        
        def step_density(x, step_height):
            """
            x is an array of positions in mm, step_height is the height of the step in nm
            """
            return (1/step_height * x / np.sqrt(crystal_curvature_radius**2-x**2)) 
        
        x = np.arange(-3, 3, 0.001) #positions in mm. The increments of 1 um are linked to the number 124 below and should not be changed
        step_density = step_density(x, step_height)
        
        som = np.cumsum(step_density)
        probed_step_density = (np.roll(som,-62)-np.roll(som,62))/124 #averages the step densities probed by the 124 um beam (so 62 um left and 62 um right of the center of the beam)
        probed_step_density = probed_step_density[63:-63] #cut off the first and last part of the array, needs to be done because the roll function does not return correct values there
        x2 = x[63:-63] 
        
        self.probed_step_density = np.interp(self.pos-center, x2, probed_step_density) #return the step density values for the measured positions
        
        self.terrace_ratio = 1/np.absolute(self.probed_step_density)/np.sqrt(3)/radius - 1
        self.terrace_fraction = self.terrace_ratio / (self.terrace_ratio+1)
        
        if plot:
            plt.plot(x2, probed_step_density)
            plt.plot(x,step_density)
            plt.scatter(self.pos-center,self.probed_step_density)
            plt.title('Step density')
            plt.xlabel('Position relative to center (mm)')
            plt.ylabel('Step density')
        #    plt.axis([-0.2,0.2,0,0.1])
            plt.show()
            plt.close()
            
        #    plt.scatter(self.pos, self.terrace_ratio)
            plt.scatter(self.pos, self.terrace_fraction*100)
            plt.title('Percentage of terraces')
            plt.xlabel('position on crystal (mm)')
            plt.ylabel('% terrace')
            plt.show()
            plt.close()
        
        if save:
            np.savetxt(savefolder+'step_density.txt',np.column_stack((self.pos,self.pos-center,self.probed_step_density)),header='position on crystal (mm), position (center at 0), step density /nm')
        
    def test_convolution(self, function='exp'):
        """
        returns the value expected for a certain function convoluted with the 
        energy distribution
        """
        if function == 'exp':
            return np.sum(self.energy_dist * np.exp(-self.energies))

        if function == 'expline':
            return np.sum(self.energy_dist * (np.exp(-self.energies)+0.05*self.energies))
        
    def linear_prediction(self, center, slopeleft, sloperight, offset, left=None, right=None):
        #print (center)
        positions = self.pos-center        
        fitleft = slopeleft * positions[left] + offset
        fitright = sloperight * positions[right] + offset
        
        return np.concatenate((fitleft,fitright))
    
    def fit_lines_without_offset(self, center,plot=False):
        """
        check what self.A_step_cross_section means exactly (I think it is a length now)
        """
        self.s0t = np.array(self.s0t)
        def fitfunction(x,a):
            return a*x
        centered_positions = self.pos-center
        left = centered_positions <= 0
        edgeleft = centered_positions > -fitrange
        left *= edgeleft
        popt, pcov = curve_fit(fitfunction, self.probed_step_density[left], self.data[left]-self.s0t[left])
        self.B_step_cross_section = popt[0]
        
        right = centered_positions > 0
        edgeright = centered_positions < fitrange
        right *= edgeright
        popt, pcov = curve_fit(fitfunction, self.probed_step_density[right], self.data[right]-self.s0t[right])
        self.A_step_cross_section = popt[0]
        
        if plot:
            plt.plot(self.probed_step_density, self.data-self.s0t, 'o',c=self.color)
            plt.plot(self.probed_step_density[left],self.probed_step_density[left]*self.B_step_cross_section, c=self.color)
            plt.plot(self.probed_step_density[right],self.probed_step_density[right]*self.A_step_cross_section, c=self.color)
        
    #calculates the cross section, and adds the probed step density for A and B to the object. Implemented after making other functions, so most functions don't use this
    def calc_cross_section(self):
        left = self.probed_step_density <= 0
        right = self.probed_step_density > 0
        self.probed_step_density_A = self.probed_step_density[right]
        self.probed_step_density_B = -self.probed_step_density[left]
        self.s0_A = self.data[right]
        self.s0_B = self.data[left]    
        self.cross_section_A = ((self.s0_A-self.s0t[right])/self.probed_step_density_A * unit_cell_width)
        self.error_cs_A = self.error[right]*unit_cell_width/self.probed_step_density_A
        self.cross_section_B = ((self.s0_B-self.s0t[left])/self.probed_step_density_B * unit_cell_width)
        self.error_cs_B = self.error[left]*unit_cell_width/self.probed_step_density_B
    
    #calculate the ratio of the A and B type cross section, average over step densities >0.18
    def calc_BAratio(self):
        datamask = self.probed_step_density_A > 0.18
        cross_section_B_interp = np.interp(self.probed_step_density_A[datamask], 
                                           np.flip(self.probed_step_density_B),
                                           np.flip(self.cross_section_B))
       # print(self.probed_step_density_A[datamask])
       # print(self.probed_step_density_B)
       # print(cross_section_B_interp)
        self.BAratio = np.average(cross_section_B_interp)/np.average(self.cross_section_A[datamask])

    def save_crossec_s0(self,folder_out):
        E = np.round(self.E_avg, 1)
        A_out = np.column_stack((self.probed_step_density_A, self.s0_A, self.cross_section_A))
        np.savetxt(folder_out+'s0_Atype_'+str(E)+'.txt', A_out, header='step density (/nm), s0, cross section (nm^2)')
        B_out = np.column_stack((self.probed_step_density_B, self.s0_B, self.cross_section_B))
        np.savetxt(folder_out+'s0_Btype_'+str(E)+'.txt', B_out, header='step density (/nm), s0, cross section (nm^2)')



    
"""
FUNCTIONS
"""  
def move_files(folder_out):
    for energy in datadic.values():
        shutil.copy(folderstart+year+datadic[energy].kw_day+file_kw, folder_out+'s0_'+str(energy)+'.txt')        



"""
FUNCTIONS FOR FIGURES
"""

#figure 4
def plot_cross_section_stepdensity(datalist,save=False,):
    #for colormap
    cmap = cm.viridis_r
    Energies = np.array([dataset.E_avg for dataset in datalist])
    max_E = np.max(Energies)
    def to_color(Energy):
        return cmap(Energy/max_E)

    #plot 
    fig, ax = plt.subplots()
    axes = [-0.65, 0.65, 0, 0.7]
    ax.axis(axes)
    ax.add_artist(Rectangle((-0.1,axes[2]),0.2,axes[3]-axes[2],color='0.95'))
    ax.text(0.8,0.9,'A type',color=cm.RdBu(0.85),transform=ax.transAxes)
    ax.text(0.1,0.9,'B type',color=cm.RdBu(0.15),transform=ax.transAxes)

    #plot data
    for dataset in datalist:
        crossec =  (dataset.data-dataset.s0t)/np.absolute(dataset.probed_step_density)*unit_cell_width
        left = dataset.probed_step_density < 0
        right = dataset.probed_step_density >= 0
        ax.plot(dataset.probed_step_density, crossec,'o-', c=to_color(dataset.E_avg))    
        # ax.plot(dataset.probed_step_density[left], crossec[left], ls='-',marker='+',markersize=7, c=to_color(dataset.E_avg))
        # ax.plot(dataset.probed_step_density[right], crossec[right], ls='-',marker='.',markersize=10, c=to_color(dataset.E_avg))

    # plt.legend(loc='upper left',bbox_to_anchor=(1,1))
    # plt.title('cross section as a function of step density')

    #make colorbar
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=0, vmax=max_E))
    cbar = fig.colorbar(sm, ax=ax)
    cbar.ax.set_ylabel('Kinetic energy ($kJ/mol$)')  

    #ticks and labels
    ticklabels = ax.get_xticks()
    ax.set_xticklabels(np.round(np.absolute(ticklabels),1)) 
    ax.set_xlabel('Step density ($nm^{-1}$)')
    ax.set_ylabel('Step reaction cross section ($nm^2$)')
    ax.tick_params(top=True, direction='in')  
    ax.tick_params(right=True, direction='in')
    
    if save:
        plt.savefig(savefolder+'cross_section_step_density.png',bbox_inches='tight',dpi=500)    
        plt.savefig(savefolder+'cross_section_step_density.eps',bbox_inches='tight',dpi=500)    
    plt.show()
    plt.close()

#for supporting information
def plot_individual_cross_section(datalist, save=False):
    ncols = 3
    nrows = int(np.ceil(len(datalist)/3))
    fig, axes2d = plt.subplots(ncols=ncols, nrows=nrows, sharex='all', sharey='row', figsize=(10,10))
    axes = axes2d.flatten()
    fig.subplots_adjust(hspace=0, wspace=0)
    for dataset, index in zip(datalist, range(len(datalist))):
        crossec =  (dataset.data-dataset.s0t)/np.absolute(dataset.probed_step_density)*unit_cell_width
        left = dataset.probed_step_density < -0.18
        right = dataset.probed_step_density >= 0.18
        # axes[index].plot(np.absolute(dataset.probed_step_density[left]), crossec[left], ls='-',marker='+',markersize=10,c=cm.RdBu(0.15)) #old way of plotting, without errorbars
        # axes[index].plot(np.absolute(dataset.probed_step_density[right]), crossec[right], ls='-',marker='.',markersize=10,c=cm.RdBu(0.85))
        axes[index].errorbar(dataset.probed_step_density_B, dataset.cross_section_B,yerr=2*dataset.error_cs_B, ls='none', capsize=5,marker='+',markersize=10,c=cm.RdBu(0.15))
        axes[index].errorbar(dataset.probed_step_density_A, dataset.cross_section_A,yerr=2*dataset.error_cs_A, ls='none', capsize=5,marker='.',markersize=10,c=cm.RdBu(0.85))
        axes[index].axis([0, 0.7, 0, 2*dataset.cross_section_B[-4]]) #not so nice way of making the y axis limits work
        axes[index].text(0.65, 0.85, str(np.round(dataset.E_avg, 1))+' kJ/mol',transform=axes[index].transAxes)

        # plt.xlabel('step density')
        # plt.ylabel('Cross section')
    for ax in axes:
        ax.tick_params(top=True, direction='in')  
        ax.tick_params(right=True, direction='in')
    
    vmiddle = int(axes2d.shape[0]/2)
    hmiddle = int(axes2d.shape[1]/2)
    axes2d[vmiddle][0].set_ylabel('Step reaction cross section ($nm^2$)')
    axes2d[axes2d.shape[0]-1][hmiddle].set_xlabel('Step density ($nm^{-1}$)')
    
    if save:
        if not os.path.exists(savefolder+'AB/'):
                os.makedirs(savefolder+'AB/')
        plt.savefig(savefolder+'AB/individual.png',bbox_inches='tight')
        plt.savefig(savefolder+'AB/individual.eps',bbox_inches='tight')
    plt.show()
    plt.close()

        

#figure 3
def plot_cross_section_E(datalist, save=True, max_E=14):
    aspectratio = 0.8
    
    #for colormap
    stepdensities = datalist[0].probed_step_density
    max_density = np.max(np.absolute(stepdensities))
    def to_color(step_density):
        return cm.RdBu(step_density/max_density*0.5+0.5)

    #plot
    fig, ax = plt.subplots()


    #insert PES
    picheight = 2*unit_cell_width #in nm, actual size of the PES
    PES = mpimg.imread(loadfolder+'PES.png')
    wpix = np.shape(PES)[1] #width of the original image in pixels
    hpix = np.shape(PES)[0]
    w = 0.65 #width of the picture relative to the figure
    h = hpix/wpix*w/aspectratio
    spacing = 0.00
    axins = ax.inset_axes([1-w-spacing*aspectratio, 1-h-spacing, w, h]) #for inset
    axins.tick_params(axis='both',bottom=False, left=False, right=False, top=False, labelbottom=False, labelleft=False)
    axins.imshow(PES)
    # ax.text(0.2,1-0.25*h-spacing*aspectratio, 'Pt(211)',transform=ax.transAxes, fontsize=9)
    ax.text(0.2,1-0.5*h-spacing*aspectratio,'Pt(211)\nTop view',transform=ax.transAxes, fontsize=9)
    
    #insert Pt(211) side view
    SV = mpimg.imread(loadfolder+'Pt211.png') 
    wpix_SV = np.shape(SV)[1] #width of the original image in pixels
    hpix_SV = np.shape(SV)[0]
    w_SV = 0.35 #width of the picture relative to the figure
    h_SV = hpix_SV/wpix_SV*w_SV/aspectratio
    spacing_SV = 0.03
    ax_SV = ax.inset_axes([1-0.5*w-spacing*aspectratio-0.5*w_SV, 1-h-spacing-h_SV-spacing_SV, w_SV, h_SV]) #for inset
    ax_SV.tick_params(axis='both',bottom=False, left=False, right=False, top=False, labelbottom=False, labelleft=False)
    ax_SV.imshow(SV)
    ax_SV.axis('off')
    ax.text(0.2,1-h-0.5*h_SV-spacing*aspectratio-spacing_SV,'Side view',transform=ax.transAxes, fontsize=9)
    dx = 0.5*(w-w_SV) 
    dy = spacing_SV + 0.25*h_SV 
    arrowlength = 0.9 #scaling factor for dx and dy
    #left arrow
    ax.arrow(1-spacing*aspectratio-w+dx,1-spacing-h-dy,-dx*arrowlength,dy*arrowlength,transform=ax.transAxes,length_includes_head=True,head_width=0.02,width=0.001,head_length=0.02,color='black')
    #right arrow
    ax.arrow(1-dx,1-spacing-h-dy,dx*arrowlength,dy*arrowlength,transform=ax.transAxes,length_includes_head=True,head_width=0.02,width=0.001,head_length=0.02,color='black')


    avg_A = [] #average cross section for A type step
    avg_B = []
    for dataset in datalist:
        left = dataset.probed_step_density <= 0
        right = dataset.probed_step_density > 0
    #    dataset.cross_section_A = ((dataset.data[right]-dataset.s0t[right])/dataset.probed_step_density[right] * unit_cell_width)
    #    dataset.cross_section_B = -((dataset.data[left]-dataset.s0t[left])/dataset.probed_step_density[left] * unit_cell_width)
    #     scatter plot
        ax.scatter(np.full(len(dataset.cross_section_B),dataset.E_avg),dataset.cross_section_B,marker='+',c=to_color(dataset.probed_step_density[left]))
        ax.scatter(np.full(len(dataset.cross_section_A),dataset.E_avg),dataset.cross_section_A,marker='.',c=to_color(dataset.probed_step_density[right]))       
    
        #calculate average for the higher energies: use stepdensity as weights because the relative error is ~1/stepdensity
        avg_A.append(np.average(dataset.cross_section_A, weights=np.absolute(dataset.probed_step_density[right])))
        avg_B.append(np.average(dataset.cross_section_B, weights=np.absolute(dataset.probed_step_density[left])))
    #plot weighted averages for each step type and energy
    
    
    #make colorbar
    sm = plt.cm.ScalarMappable(cmap=cm.RdBu, norm=plt.Normalize(vmin=-max_density, vmax=max_density))
    cbar = fig.colorbar(sm, ax=ax)
    labels = cbar.ax.get_yticks()

    #axes and labels
    ax.tick_params(top=True, direction='in')  
    ax.tick_params(right=True, direction='in')
    labels = cbar.ax.get_yticks()
    cbar.ax.set_yticklabels(np.round(np.absolute(labels),1))
    cbar.ax.set_ylabel('step density ($nm^{-1}$)')   
    plt.axis([0,max_E,0,0.6])

    ax.set_aspect(aspect=aspectratio/ax.get_data_ratio())
    plt.xlabel('Kinetic energy ($kJ/mol$)')
    plt.ylabel('Step reaction cross section ($nm^2$)')

   


    #plot horizontal lines
    avg_highsd = np.average(np.array([x.cross_section_A[-1] for x in datalist])[7:])
    avg_lowsd = np.average(np.array([x.cross_section_A[-8] for x in datalist])[7:])
    print('lowsd ',datalist[0].probed_step_density[-8])

    col_indirect = [0.35,0.35,0.35]

    ax.plot([0, max_E], [avg_lowsd, avg_lowsd], '--', color='black')
    ax.plot([0, max_E], [avg_highsd, avg_highsd], '-', color='black')
    ax.plot([0, datalist[0].E_avg], [datalist[0].cross_section_A[-8],datalist[0].cross_section_A[-8]],'--',color=col_indirect)
    ax.plot([0, datalist[0].E_avg], [datalist[0].cross_section_A[-1],datalist[0].cross_section_A[-1]],color=col_indirect)
    print('cross section lowsd ',avg_lowsd, ', cross section highsd ', avg_highsd,', atom area ',np.pi*(unit_cell_width/2)**2)

    #insert circles
    r = np.sqrt(avg_lowsd/np.pi) / picheight * hpix #radius of the circle in the inset
    circ = Circle((0.51*wpix, 0.5*hpix), radius=r, fill=False, ls='--', color='black')
    axins.add_artist(circ)
    r = np.sqrt(avg_highsd/np.pi) / picheight * hpix #radius of the circle in the inset
    circ = Circle((0.51*wpix, 0.5*hpix), radius=r, fill=False, color='black')
    axins.add_artist(circ)

    #insert rectangles
    h = unit_cell_width / picheight * hpix
    w = (datalist[0].cross_section_A[-1]-avg_highsd) / unit_cell_width / picheight * hpix
    print('area indirect mechanism ',datalist[0].cross_section_A[-1]-avg_highsd)
    axins.add_artist(Rectangle((wpix*0.65-0.5*w+r,0.25*hpix),w,h,fill=False,color=col_indirect))

    w = (datalist[0].cross_section_A[-8]-avg_lowsd) / unit_cell_width / picheight * hpix
    print('area indirect mechanism ',datalist[0].cross_section_A[-8]-avg_lowsd)
    axins.add_artist(Rectangle((wpix*0.65-0.5*w+r,0.25*hpix),w,h,fill=False,color=col_indirect,ls='--'))


    
    if save:
        plt.savefig(savefolder+'cross_section.png',dpi=500)
        plt.savefig(savefolder+'cross_section.eps',dpi=500, bbox_inches='tight')
    
    plt.show()
    plt.close()
    

def plot_s0_vs_pos(datalist, save=False):
    #plot all s0 in a single plot
    for dataset in datalist:
        dataset.plot_s0()
    plt.axis([17.9,21.9,0,None])
    plt.xlabel('Position (mm)')
    plt.ylabel('s0')
    plt.legend(loc='upper left',bbox_to_anchor=(1,1))
    if save:
        plt.savefig(savefolder+'all_s0.png',dpi=500,bbox_inches='tight')
    plt.show()
    plt.close()     

#figure 2
def plot_s0_vs_E(datalist, use_stepdens=True, save=False, max_E=14): 
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
        dataset.s0_const = [] #empty list for constant contribution to s0
        dataset.popt = [] #empty list for all fit parameters
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
        return cm.RdBu(step_density/max_density*0.5+0.5) #*0.5 + 0.5 to take only one half of the colormap for positive step densities, and one half for negative step densities
    
    #plot
    fig, axes = plt.subplots(4,1,sharex='col',figsize=(8,10))
    ax1 = axes[0]
    ax2 = axes[1]
    ax3 = axes[2]
    ax4 = axes[3]
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
    
    #make array of positions/step densities to iterate over, make sure it is in the right order for plotting
    # positions = np.array(list(dic.keys()))
    # sort = np.argsort(np.absolute(positions))
    # positions = np.flip(positions[sort])
    for k in dic.keys():
        #FIT
        popt, pcov = curve_fit(expline, E, dic[k],p0=guess, bounds=(lower,upper))
        #ADD S0T TO DATA OBJECTS
        for dataset in datalist:
            dataset.s0t.append(popt[2]*dataset.E_avg)
            dataset.s0_const.append(popt[3])
            dataset.popt.append(popt)
        
        #PLOT ORIGINAL AND FIT
        if k<=0:
            # ax2.plot(E,popt[2]*E+popt[3],color=to_color(k)) #plot just the line
            ax2.plot(E, dic[k],'o', color=to_color(k)) 
            ax2.plot(Efit, expline(Efit,*popt),'-',label=str(k), color=to_color(k))
        else:
            ax1.plot(E, dic[k],'o',label=str(k), color=to_color(k))
            ax1.plot(Efit, expline(Efit,*popt),'-',label=str(k), color=to_color(k))

    #make np array instead of list 
    for dataset in datalist:
        dataset.s0t=np.array(dataset.s0t)   
        dataset.s0_const=np.array(dataset.s0_const)

    #plot original and subtracted data
    pos_index = 0
    for dataset in datalist:
        ax3.plot(dataset.E_avg, dataset.data[pos_index],'o', color='black')
        ax3.plot(dataset.E_avg, dataset.data[pos_index]-dataset.s0t[pos_index], 'o', color='gray')
    ax3.plot(E, expline(E,*dataset.popt[pos_index]), '--', c='black')
    ax3.plot(E, dataset.s0_const[pos_index]*np.ones(len(E)),color='gray')
    ax3.plot(E, dataset.s0t[pos_index]/dataset.E_avg*E+dataset.s0_const[pos_index], color='black')
    #add an arrow
    arr = Arrow(12,0.2,0,-0.1,color='black',width=0.2)
    ax3.add_artist(arr)


    #plot the energy distributions
    for dataset in datalist:
        ax4.plot(dataset.energies, dataset.energy_dist, color=cm.gray(dataset.E_avg/max_E))
    # aspectratio = 0.05
    # ax4.set_aspect(aspect=aspectratio/ax4.get_data_ratio())
    
    #axis limits
    axlims = [0,max_E,0,0.45]
    for ax in axes[:3]:
        ax.axis(axlims)
    ax4.axis([0,max_E,0,None])   
    
    #make colorbar
    sm = plt.cm.ScalarMappable(cmap=cm.RdBu, norm=plt.Normalize(vmin=-max_density, vmax=max_density))
    #cbar = fig.colorbar(sm, ax=[ax1, ax2])
    cbar = fig.colorbar(sm, ax=axes.ravel().tolist(),shrink=0.5,anchor=(0,1))

    labels = cbar.ax.get_yticks()
    cbar.ax.set_yticklabels(np.round(np.absolute(labels),1))
    cbar.ax.set_ylabel('step density ($nm^{-1}$)') 


    #ticks
    for ax in axes:
        ax.tick_params(top=True, direction='in')  
        ax.tick_params(right=True, direction='in')
    ax4.tick_params(axis='y', left=False, labelleft=False, right=False)

    #labels and title
    for ax,letter,text in zip(axes,'abcd',['A type', 'B type', '','']):
        ax.text(0.012,0.9,letter,transform=ax.transAxes, fontsize=14)
        ax.text(0.83,0.87,text,transform=ax.transAxes, fontsize=16)        
    ax4.set_xlabel('Energy ($kJ/mol$)')
    for ax in axes[:3]:
        ax.set_ylabel('$s_0$')
    ax4.set_ylabel('Beam distribution')
    # ax1.set_title('B type (red) and A type (blue)')
    if save:
        plt.savefig(savefolder+'s0_vs_E.png',dpi=500)
        plt.savefig(savefolder+'s0_vs_E.eps',dpi=500, bbox_inches='tight')
    plt.show()        
 
#supporting information    
def test_convolution(datalist, save=False, max_E=14):
    fig, axes = plt.subplots(nrows=2, sharex=True)
    fig.subplots_adjust(hspace=0)
    #plot exp+line
    x = np.arange(0.7,13.1, 0.01)
    axes[0].plot(x, np.exp(-x)+0.05*x,c='black', label='Function without convolution')
    #test convolution
    for dataset in datalist:
        axes[0].plot(dataset.E_avg, dataset.test_convolution(function='expline'), 'o',c=str(dataset.E_avg/max_E), label='Convolution')


    for dataset in datalist:
        axes[1].plot(dataset.energies, dataset.energy_dist, color=str(dataset.E_avg/max_E))
    
    #labels etc
    for ax in axes:
        ax.tick_params(axis='y', left=False, right=False, labelleft=False)
        ax.axis([0,max_E,0,None])
    axes[0].set_ylabel('s0')
    axes[1].set_ylabel('Beam distribution')
    axes[1].set_xlabel('Energy ($kJ/mol$)')
    # axes[0].legend()


    
    if save:
        plt.savefig(savefolder+'convolution.png',dpi=500)
        plt.savefig(savefolder+'convolution.eps',dpi=500)
    plt.show()
    plt.close()   

#supporting information
def plot_vs_temperature(save=False):
    path1 = folderstart+"2020/07 Jul/200707/KW/Images/s0_pos_temp.txt"
    path2 = folderstart+"2020/07 Jul/200709/KW/Images/s0_pos_temp.txt"
    data1 = np.loadtxt(path1, skiprows=1,unpack=False)
    data2 = np.loadtxt(path2, skiprows=2,unpack=False)
    data = np.concatenate((data1,data2), axis=0)
    temp = data[:,0]
    pos = data[:,1]
    error = data[:,3]
    data = data[:,2]

    #convert positions to step densities
    pos_to_stepdens={18:0.59,19:0.29,21:0.3,22:0.6}

    #for colors
    def to_color(pos):
        return cm.RdBu((pos-17.5)/5)
    def marker(pos):
        if pos < 20: return '+'
        if pos >= 20: return '.'
    for position in np.unique(pos):
        index = pos == position
        plt.errorbar(temp[index],data[index],yerr=2*error[index],capsize=4,ls='None',marker=marker(position),markersize=10,label=str(pos_to_stepdens[position]),c=to_color(position))
    plt.legend()
    plt.axis([90,260,0,0.5])
    plt.xlabel('Temperature ($K$)')
    plt.ylabel('$s_0$')
    if save:
        plt.savefig(savefolder+'temperature.png',dpi=500)  
        plt.savefig(savefolder+'temperature.eps',dpi=500)  

    plt.show()
    plt.close()


#supporting? comparison to Irenes data
def compare_irene(datalist,save=False):
    stepdensity_irene = [0.73,1.13,1.56]
    s0_irene = [0.12, 0.14, 0.205]
    plt.plot(stepdensity_irene, s0_irene,'*',c=cm.RdBu(0.85),label='Groot et al., 2013 (A type)')
    A = datalist[0].probed_step_density > 0 
    B = datalist[0].probed_step_density <= 0 
    plt.plot(np.absolute(datalist[0].probed_step_density[B]),datalist[0].s0_const[B],'+',c=cm.RdBu(0.15),label='B type')
    plt.plot(datalist[0].probed_step_density[A],datalist[0].s0_const[A],'.',c=cm.RdBu(0.85),label='A type')
    plt.axis([0, 1.6,0,0.25])
    plt.xlabel('step density ($nm^{-1}$)')
    plt.ylabel('s0 intercept')
    plt.legend()
    if save:
        plt.savefig(savefolder+'intercept.png',dpi=500)  
        plt.savefig(savefolder+'intercept.eps',dpi=500)         
    plt.show()
    plt.close()


#supporting, fig 3 + hard cube
def plot_cross_section_E_hardcube(datalist,Ewell_list, save=True, max_E=14):
    
    #for colormap
    stepdensities = datalist[0].probed_step_density
    max_density = np.max(np.absolute(stepdensities))
    def to_color(step_density):
        return cm.RdBu(step_density/max_density*0.5+0.5)

    #plot
    fig, axes = plt.subplots(len(Ewell_list),sharex=True)
    fig.subplots_adjust(hspace=0)

    #make single array with cross sections
    energies_1D = [dataset.E_avg for dataset in datalist]
    stepdens_1D = datalist[0].probed_step_density
    energies, stepdens = np.meshgrid(energies_1D, stepdens_1D)
    crossecs = np.zeros(np.shape(energies))
    for i in range(len(datalist)):
        crossecs[:,i] = (datalist[i].data-datalist[i].s0t)/datalist[i].probed_step_density*unit_cell_width
 

    for ax,Ewell in zip(axes,Ewell_list):
        #load hard cube data
        E_hardcube, trap_hardcube_abs = np.loadtxt(loadfolder+'hardcube_'+str(Ewell)+'.txt',skiprows=1,unpack=True)

        # for i in range(np.shape(crossecs)[1]):
        for i in [0, 5, -6, -1]:
            if i < 0:
                marker = '.'
            else:
                marker= '+'
            ax.plot(energies[i,:],np.absolute(crossecs[i,:]),marker,c=to_color(stepdens[i][0]))
            #scale hard cube 
            avg_crossec = np.average(np.absolute(crossecs[i,7:]))
            trap_hardcube = trap_hardcube_abs / np.max(trap_hardcube_abs) * (np.absolute(crossecs[i,0]) - avg_crossec) + avg_crossec
            ax.plot(energies[i,:],trap_hardcube,linewidth=1,c=to_color(stepdens[i][0]))
        
        #axes and labels
        ax.text(0.6,0.8,str(Ewell)+' kJ/mol well depth',transform=ax.transAxes)
        ax.tick_params(top=True, direction='in')  
        ax.tick_params(right=True, direction='in')
        ax.axis([0,max_E,0,0.25])
    
    ax.set_xlabel('Kinetic energy ($kJ/mol$)')

    fig.add_subplot(111,frameon=False)
    plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
    plt.ylabel('Step reaction cross section ($nm^2$)')


    #make colorbar
    sm = plt.cm.ScalarMappable(cmap=cm.RdBu, norm=plt.Normalize(vmin=-max_density, vmax=max_density))
    cbar = fig.colorbar(sm, ax=axes.ravel().tolist())
    labels = cbar.ax.get_yticks()
    labels = cbar.ax.get_yticks()

    cbar.ax.set_yticklabels(np.round(np.absolute(labels),1))
    cbar.ax.set_ylabel('step density ($nm^{-1}$)')   


    if save:
        plt.savefig(savefolder+'cross_section_hardcube.png',dpi=500)
        plt.savefig(savefolder+'cross_section_hardcube.eps',dpi=500, bbox_inches='tight')
    
    plt.show()
    plt.close()

#calculate area of attractive potential
def calc_attr():
    PES = mpimg.imread(loadfolder+'PES2_area.png')
    wpix = np.shape(PES)[1]
    hpix = np.shape(PES)[0]
    print(np.shape(PES))
    print(PES[150,500])
    flattened = np.average(PES,axis=2)
    print(np.shape(flattened))
    a = np.argwhere(np.absolute(flattened-0.25) < 0.05)
    print(np.shape(a))
    area_fraction = np.shape(a)[0]/wpix/hpix*4
    print('area_fraction in PES',area_fraction)
    area_nm2 = area_fraction*unit_cell_width*unit_cell_width*wpix/hpix
    print('area in nm2 in PES', area_nm2)
    print('area in nm2 in data', 0.15)
    print('area fraction in data', 0.15/unit_cell_width*0.57)#0.57 is the step density



def plot_energy_dist(datalist):
    #plot all energy distributions in a single plot
    for dataset in datalist:
        dataset.plot_TOF()
    plt.axis([0,15,0,None])
    plt.xlabel('Energy (kJ/mol)')
    plt.ylabel('Probability density')
    plt.title('Energy distributions of measured beams')
    plt.legend()
    plt.savefig(savefolder+'energy_distributions.png')
    plt.show()
    plt.close()



"""
FIT FUNCTION
"""
def fit_all_slopes(datalist, fitrange=2.0, center_guess=20, no_offset=False,save=False, plot=False):
    def residuals(fitparams):
        """
        fitparams = [center, slopeleft1, sloperight1, offset1, slopeleft2, sloperight2, etc]
        fits as if the center is at x=0
        use_corrected: instead of raw data, fit only step contribution
        """
        res = np.array([])
        concatenated_data=np.array([])
       # print(fitparams[0])
        for dataset,i in zip(datalist,range(0,3*len(datalist),3)):
            
            positions = dataset.pos-center_guess
            left = positions<0
            right = np.invert(left)
            edgeleft = positions > -fitrange
            edgeright = positions < fitrange
            left *= edgeleft
            right *= edgeright
            concatenated_data = np.concatenate((concatenated_data, dataset.data[left+right]))

            res = np.concatenate((res, dataset.linear_prediction(fitparams[0],
                                        *fitparams[i+1:i+4], left=left, right=right)))
        assert len(res) == len(concatenated_data)
        return res - concatenated_data   

    def guess_fitparams(center_guess):
        fitparams=np.array([center_guess])
        for dataset in datalist:
            positions = dataset.pos-center_guess
            left = positions<0
            right = np.invert(left)
            edgeleft = positions > -fitrange
            edgeright = positions < fitrange
            left *= edgeleft
            right *= edgeright

            slopeleft, offsetleft = np.polyfit(positions[left],dataset.data[left],1)
            sloperight, offsetright = np.polyfit(positions[right],dataset.data[right],1)            
            offset = (offsetleft+offsetright)/2
            fitparams = np.concatenate((fitparams,np.array([slopeleft, sloperight, offset])))

        return fitparams
    
    concatenated_data = np.array([])
    for dataset in datalist:
        concatenated_data = np.concatenate((concatenated_data, dataset.data))
    
    fitparams = guess_fitparams(center_guess)
   # print (fitparams,'\n')
    upper = fitparams+np.absolute(fitparams)*0.1
    lower = fitparams-np.absolute(fitparams)*0.1
    if no_offset:
        fitparams[3::3]=0
        upper[3::3]=1E-5
        lower[3::3]=0

    bounds = (lower, upper)
    
    
    result = least_squares(residuals, fitparams, bounds=bounds)
   # print(result.x)
    print('center of crystal at '+str(result.x[0])+' mm')

    #Plot the fitted lines in the sticking data, and save fitted lines
    for dataset,i in zip(datalist, range(0,3*len(datalist),3)):
        center = result.x[0] 

        positions = dataset.pos-center #convert positions so center is at x=0
        left = positions<0 #data mask for points left of center
        right = np.invert(left)
        edgeleft = positions > -fitrange #data mask for points inside fit range
        edgeright = positions < fitrange
        left *= edgeleft #combine two data masks to create final mask
        right *= edgeright
        posleft = dataset.pos[left]
        posright = dataset.pos[right]
        dataleft = result.x[i+1]*positions[left]+result.x[i+3]
        dataright = result.x[i+2]*positions[right]+result.x[i+3]
        if save:
            np.savetxt(savefolder+str(np.round(dataset.E_avg,decimals=1))+'_fitleft_range_'
                       +str(fitrange)+'.txt', np.column_stack((posleft-center, dataleft)), header='position left (mm), s0 fit left')
            np.savetxt(savefolder+str(np.round(dataset.E_avg,decimals=1))+'_fitright_range_'
                       +str(fitrange)+'.txt', np.column_stack((posright-center, dataright)), header='position right (mm), s0 fit right')
        
        if plot:
            plt.plot(posleft, dataleft, c=dataset.color)
            plt.plot(posright, dataright, c=dataset.color)    
            plt.scatter(dataset.pos,dataset.data, c=dataset.color)
    if plot:
        plt.axis([None,None,0,0.5])
        plt.savefig(savefolder+'fitted.png',dpi=500)
        plt.show()
        plt.close()
    
    center=result.x[0]
    slopeleft=result.x[1::3]
    sloperight=result.x[2::3]
    offset=result.x[3::3]

    return center, slopeleft, sloperight, offset


     
main()