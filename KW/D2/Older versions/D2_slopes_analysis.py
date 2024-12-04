# -*- coding: utf-8 -*-
"""
Created on Fri Feb 28 15:05:09 2020

@author: jansenc3
"""

import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import least_squares
plt.style.reload_library()
#plt.style.use('voorbeeld')

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
           11.0:{'kw_day':'07 Jul/200727', 'TOF_day':'07 Jul/200727','c':'black'}}

datalist=[] #will be filled with data objects


fitrange = 1.2 #program fails for some values due to arrays not being the same size

def main():
    for dataset in datadic.values():
        datalist.append(data(folderstart, year, dataset['kw_day'],dataset['TOF_day'], 
                             file_kw, file_TOF, dataset['c'])) #list of data objects
     

    
    plot_energy_test_convolution(datalist)
    
    #plot all s0 in a single plot
    for dataset in datalist:
        dataset.plot_s0()
    plt.axis([17.9,21.9,0,None])
    plt.legend()
    plt.show()
    plt.close()     
    
    
    center = fit_all_slopes(datalist,fitrange=fitrange,center_guess=19.9)
    
    dataset.convert_step_density(center, save=True)

    
    
    
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
    
        

    def convert_step_density(self, center,save=False):
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
        




def fit_all_slopes(datalist, fitrange=2.0, center_guess=20):
    def residuals(fitparams):
        """
        fitparams = [center, slopeleft1, sloperight1, offset1, slopeleft2, sloperight2, etc]
        fits as if the center is at x=0
        fitrange should be 1.0, 1.2 etc (not 1.1 etc) for some reason
        """
        res = np.array([])
        concatenated_data=np.array([])
#        print(fitparams[0])
        for dataset,i in zip(datalist,range(0,3*len(datalist),3)):
            
            positions = dataset.pos-fitparams[0]
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
    np.savetxt(savefolder+'slope_and_111_range_'+str(fitrange)+'.txt', np.column_stack((E,-slopeleft,sloperight, offset)), 
               header='Energy, B slope, A slope, s0 (111)')
    plt.plot(E,-slopeleft,label='B type', c='r',marker='o')
    plt.plot(E,sloperight,label='A type', c='b',marker='o')
    plt.plot(E,offset,label='S0 (111)', c='black',marker='o')  
    plt.legend()
    plt.savefig(savefolder+'slopes_and_s0.png',dpi=500)
    plt.show()
    plt.close()
    return center

        
        
def fitfunction(x, E, fitparams, ):
    
    def fraction_terraces(x):  # maybe better outside fitfunction, part of data
        pass
    
    def fslope(E, fitparams):
        pass
    
    def f111(E, fitparams):
        return fitparams[0]  # is maar een voorbeeld

    return fslope(E, fitparams) * x + f111(E, fitparams) * fraction_terraces(x)
        
  
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
    
        
main()