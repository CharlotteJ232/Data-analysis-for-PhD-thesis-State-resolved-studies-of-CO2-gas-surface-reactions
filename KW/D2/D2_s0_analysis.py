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
plt.style.reload_library()
#plt.style.use('voorbeeld')

unit_cell_width = 2*0.13874 #Pt atom diameter

folderstart = 'P:/Surfdrive/DATA/'
folderstart = 'D:/Surfdrive/DATA/'
year = '2020/'
file_kw = '/KW/Images/stickingprob_vs_position.txt'
file_TOF = '/TOF/Images/4.0/energy.txt'
savefolder = folderstart+'D2 sticking probability as a function of step density and kinetic energy/Images/'


datadic = {0.7:{'kw_day':'03 Mar/200310', 'TOF_day':'03 Mar/200310','c':'red'},
           0.8:{'kw_day':'05 May/200512', 'TOF_day':'01 Jan/200128','c':'orange'},
           1.2:{'kw_day':'03 Mar/200305', 'TOF_day':'02 Feb/200220','c':'yellow'},
           2.0:{'kw_day':'06 Jun/200602', 'TOF_day':'05 May/200511','c':'green'},
           2.9:{'kw_day':'03 Mar/200303', 'TOF_day':'02 Feb/200218','c':'blue'},
           4.0:{'kw_day':'03 Mar/200311', 'TOF_day':'03 Mar/200306','c':'indigo'},
           4.77:{'kw_day':'05 May/200519', 'TOF_day':'05 May/200519','c':'purple'},
           6.1:{'kw_day':'05 May/200520', 'TOF_day':'02 Feb/200204','c':'violet'},
           7.8:{'kw_day':'08 Aug/200817', 'TOF_day':'08 Aug/200813_2','c':'black'},
           9.4:{'kw_day':'08 Aug/200814', 'TOF_day':'08 Aug/200813','c':'grey'},
           10.7:{'kw_day':'07 Jul/200727', 'TOF_day':'07 Jul/200727','c':'cyan'},
           12.9:{'kw_day':'09 Sep/200904', 'TOF_day':'09 Sep/200904','c':'magenta'}}
           # 11.2:{'kw_day':'08 Aug/200825', 'TOF_day':'08 Aug/200825','c':'magenta'}}

datalist=[] #will be filled with data objects


fitrange = 0.8 #program fails for some values due to arrays not being the same size

def main():
    #read data and make data objects
    for dataset in datadic.values():
        datalist.append(data(folderstart, year, dataset['kw_day'],dataset['TOF_day'], 
                             file_kw, file_TOF, dataset['c'])) #list of data objects
    
#    look at the energy convolution in the data. May be useful to repeat this for an exponential + line later
    plot_energy_test_convolution(datalist)
    
    #plot the raw s0 data
    plot_s0_vs_pos(datalist,save=True)

    #Fit straight lines through the data, now only mostly useful to fit the center
    center,slopeleft,sloperight,offset = fit_all_slopes(datalist,fitrange=fitrange,center_guess=19.91,save=False)
    
    #make the dataset.probed_step_density array
    for dataset in datalist:
        dataset.convert_step_density(center, plot=False, save=False)
    
    #plot s0 vs E and fit the exp + line 
    plot_s0_vs_E(datalist, use_stepdens=True, save=True)   

   
    #plot subtracted dataset
    for dataset in datalist:
        plt.plot(dataset.pos, dataset.data-dataset.s0t, 'o-', c=dataset.color)
    plt.axis([None,None,0,None])
    plt.show()
    plt.close()
    
    #plot cross section for each position individually
    for dataset in datalist:
        left = dataset.probed_step_density <= 0
        right = dataset.probed_step_density > 0
        cross_section_A = ((dataset.data[right]-dataset.s0t[right])/dataset.probed_step_density[right] * unit_cell_width)
        cross_section_B = -((dataset.data[left]-dataset.s0t[left])/dataset.probed_step_density[left] * unit_cell_width)
        # remove highest and lowest data point
#        cross_section_A = np.delete(cross_section_A, [np.argmin(cross_section_A),np.argmax(cross_section_A)])
#        cross_section_B = np.delete(cross_section_B, [np.argmin(cross_section_B),np.argmax(cross_section_B)])
        #scatter plot
        plt.scatter(np.full(len(cross_section_B),dataset.E_avg),cross_section_B,marker='+',c=[((x+1)/len(cross_section_B),0,0) for x in range(len(cross_section_B))])
        plt.scatter(np.full(len(cross_section_A),dataset.E_avg),cross_section_A,marker='+',c=[(0,0,1-x/len(cross_section_A)) for x in range(len(cross_section_A))])       
        # plt.errorbar(dataset.E_avg, np.average(cross_section),yerr=2*np.std(cross_section),fmt='o',capsize=5)
        # plt.errorbar(dataset.E_avg, np.average(cross_section_A),yerr=[[np.min(cross_section_A)],[np.max(cross_section_A)]],fmt='o',capsize=5,color='blue')
        # plt.errorbar(dataset.E_avg, np.average(cross_section_B),yerr=[[np.min(cross_section_B)],[np.max(cross_section_B)]],fmt='o',capsize=5,color='red')
        
        # plt.boxplot(cross_section,positions=[dataset.E_avg])
    plt.axis([0,12,-0.25,0.65])
    plt.xlabel('Kinetic energy (kJ/mol)')
    plt.ylabel('Step reaction cross section')
    plt.show()
    plt.close()
    
    print('approx active site area radius ',np.sqrt(0.05/np.pi), ', atom radius ',unit_cell_width/2)
    
    #plot cross section calculated by fitting lines
    A_type = []
    B_type = []
    E = []
    for dataset in datalist:
        #plot the subtracted s0 data with fitted lines first
        dataset.fit_lines_without_offset(center,plot=True)
        E.append(dataset.E_avg)
        A_type.append(dataset.A_step_cross_section)
        B_type.append(dataset.B_step_cross_section)
    plt.show()
    plt.close()
    
    A_cross_section = np.array(A_type)*unit_cell_width
    B_cross_section = np.absolute(B_type)*unit_cell_width
    
    savearray = np.column_stack((E, A_cross_section,B_cross_section))
    np.savetxt(savefolder+'reaction_cross_section.txt',savearray, header='E (kJ/mol), A (nm^2), B (nm^2)')
    
    #plot the steps cross section as a function of energy 
    plt.plot(E,A_cross_section,'o-',c='b')
    plt.plot(E,B_cross_section,'o-',c='r')
    plt.axis([0,12,0,0.3])    
    plt.xlabel('Energy (kJ/mol)')
    plt.ylabel('Reaction cross section (nm^2)')
    plt.savefig(savefolder+'cross_section.png',dpi=500)
    plt.show()
    plt.close()    
        
    

    
    
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
        radius = 0.13874 #nm, for platinum
        
        def step_density(x, step_height):
            return (1/step_height * x / np.sqrt(225-x**2)) #np.absolute(x)
        
        x = np.arange(-3, 3, 0.001)
        step_density = step_density(x, step_height)
        
        som = np.cumsum(step_density)
        
        probed_step_density = (np.roll(som,-62)-np.roll(som,62))/124
        probed_step_density = probed_step_density[63:-63]
        x2 = x[63:-63]
        
        self.probed_step_density = np.interp(self.pos-center, x2, probed_step_density)
        
        self.terrace_ratio = 1/np.absolute(self.probed_step_density)/np.sqrt(3)/radius - 1
        self.terrace_fraction = self.terrace_ratio / (self.terrace_ratio+1)
        
        if plot:
            plt.plot(x2, probed_step_density)
            plt.plot(x,step_density)
            plt.scatter(self.pos-center,self.probed_step_density)
            plt.title('Step density')
            plt.xlabel('Position relative to center (mm)')
            plt.ylabel('Step density')
    #        plt.axis([-0.2,0.2,0,0.1])
            plt.show()
            plt.close()
            
    #        plt.scatter(self.pos, self.terrace_ratio)
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
        
    def linear_prediction(self, center, slopeleft, sloperight, offset, left=None, right=None):
#        print (center)
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
        
    
"""
FUNCTIONS
"""

    
"""
PLOT FUNCTIONS 
"""
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
 
    
def plot_energy_test_convolution(datalist):
    #test convolution
    for dataset in datalist:
        plt.scatter(dataset.E_avg, dataset.test_convolution(), c=dataset.color, label='Convolution')
    plt.plot(np.arange(0.7,11.1,0.01), np.exp(-np.arange(0.7,11.1,0.01)),c='black', label='Function without convolution')
    plt.title('Effect of convolution in exponential function')
    plt.xlabel('Energy (kJ/mol)')
    plt.legend()
    plt.savefig(savefolder+'convolution.png')
    plt.show()
    plt.close()   
    
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
FIT FUNCTIONs
"""
def fit_all_slopes(datalist, fitrange=2.0, center_guess=20, no_offset=False,save=False):
    def residuals(fitparams):
        """
        fitparams = [center, slopeleft1, sloperight1, offset1, slopeleft2, sloperight2, etc]
        fits as if the center is at x=0
        use_corrected: instead of raw data, fit only step contribution
        """
        res = np.array([])
        concatenated_data=np.array([])
#        print(fitparams[0])
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
#    print (fitparams,'\n')
    upper = fitparams+np.absolute(fitparams)*0.1
    lower = fitparams-np.absolute(fitparams)*0.1
    if no_offset:
        fitparams[3::3]=0
        upper[3::3]=1E-5
        lower[3::3]=0

    bounds = (lower, upper)
    
    
    result = least_squares(residuals, fitparams, bounds=bounds)
#    print(result.x)
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
        plt.plot(posleft, dataleft, c=dataset.color)
        plt.plot(posright, dataright, c=dataset.color)    
        plt.scatter(dataset.pos,dataset.data, c=dataset.color)
    plt.axis([None,None,0,0.5])
    plt.savefig(savefolder+'fitted.png',dpi=500)
    plt.show()
    plt.close()
    
    center=result.x[0]
    slopeleft=result.x[1::3]
    sloperight=result.x[2::3]
    offset=result.x[3::3]
    
    #Plot the slopes and s0(111) and save
    E=list(datadic.keys())
    if save:
        np.savetxt(savefolder+'slope_and_111_range_'+str(fitrange)+'.txt', np.column_stack((E,-slopeleft,sloperight, offset)), 
                   header='Energy, B slope, A slope, s0 (111)')
    plt.plot(E,-slopeleft,label='B type', c='r',marker='o')
    plt.plot(E,sloperight,label='A type', c='b',marker='o')
    plt.plot(E,offset,label='S0 (111)', c='black',marker='o')  
    plt.legend()
    plt.savefig(savefolder+'slopes_and_s0.png',dpi=500)
    plt.show()
    plt.close()
    return center, slopeleft, sloperight, offset




    
        
main()