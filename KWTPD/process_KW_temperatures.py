# -*- coding: utf-8 -*-
"""
Created on Mon May 20 16:30:21 2019

@author: jansenc3
"""

import numpy as np
from matplotlib import pyplot as plt
import os

folderstart = 'P:/Surfdrive/DATA/'
#folderstart = 'D:/Surfdrive/DATA/'

positions_odd = np.array([18,19,21,22])
positions_even = np.array([22,21,19,18])
positions = positions_odd
# positions_odd = np.arange(179,220, 2)
# positions_even = np.arange(219, 178, -2)
# positions=positions_odd


folder = folderstart+'2020/07 Jul/200707/KW/'
datasets = np.arange(3,20,1)
temperatures = [100,150,200,250,300,350]

folder = folderstart+'2020/07 Jul/200709/KW/'
datasets = np.arange(3,20,1)
temperatures = [100,250]

folder = folderstart+'2020/07 Jul/200713/KW/'
datasets = np.arange(1,12.1,1)
temperatures = [100,250]
positions_odd = np.array([16.5,17.5,19,21,22.5,23.5])
positions_even = np.array([23.5,22.5,21,19,17.5,16.5])
positions = positions_odd




filenamestart = 'KW' #'KW' for standard 'KW01' filenames
ext = '.txt'
savefolder = folder+'Images/'

t_flag1_open = 9 #s, total time flag 1 is open
t_open_flag2 = 3 #s, time after opening of flag 1
#t_dose = 3 #s, dose time of the beam 
t_background = 5 #amount of seconds left and right of the measurement that are kept in the datasets. this can be chosen


save_individual = False



def main():
    measurements = {}
    for temp in temperatures:
        measurements[temp]={}
        for pos in positions:
            measurements[temp][pos] = measurement(pos, t_background, t_flag1_open, t_open_flag2)
    
        read_data_files(measurements[temp], numbers='odd',filenamestart=str(temp)+'_KW') #only reads the first column of data now! add: filenamestart=str(temp)+'_KW'
        read_data_files(measurements[temp], numbers='even',filenamestart=str(temp)+'_KW')
        
        
    stickingdic = {}
    for position in positions:
        stickingdic[position] = []
    
    for temp in temperatures:
    # errorlist = []
        for position in positions:
            measurements[temp][position].plot_measurement_set(save=save_individual)    
            measurements[temp][position].sum_measurements()
            measurements[temp][position].flip_normalize()
            measurements[temp][position].fit_stickingprob(t_fit=2.2,t_shift=0.5)
            measurements[temp][position].plot_analyzed_data(save=save_individual)
            stickingdic[position].append(measurements[temp][position].stickingprob)
            # errorlist.append(measurements[position].sigma)
    
    for position in positions:

        
        for temp in temperatures:
            plt.plot(list(measurements[temp][position].timedict.values())[0],
                     measurements[temp][position].normalizeddata, label=str(temp)+'K')
        plt.title(str(position)+' mm')
        plt.axis([7,13,-0.1,0.7])
        plt.xlabel('Time (s)')
        plt.ylabel('Sticking probability')
        plt.legend()
        plt.savefig(savefolder+str(position)+'.png',dpi=500)
        plt.show()

    for position in positions:
        center = 20
        col = 1 - np.absolute(position-center) / np.max(positions-center)
        if position-center<0:
            color = (1, col, col)
        else:
            color = (col, col, 1)
            
        stickingdic[position] = np.array(stickingdic[position])
        plt.plot(temperatures, stickingdic[position],'o-',c=color,label=str(position)+' mm')
    plt.xlabel('Temperature(K)')
    plt.ylabel('S0')
    plt.legend()
    plt.savefig(savefolder+'all.png',dpi=500)
    plt.show()
    plt.close()

########## End main ######################
    
########## Classes  ######################
    
class measurement:
    def __init__(self, position, t_background, t_flag1_open, t_open_flag2):
        self.position = position
        self.t_background = t_background
        self.t_flag1_open = t_flag1_open
        self.t_open_flag2 = t_open_flag2

        self.timedict = {}
        self.datadict = {}
  
    
    def plot_measurement_set(self, save=False):
        for name in self.datadict.keys():
            measurement_number = float(name.replace(filenamestart,''))         
            plt.plot(self.timedict[name],self.datadict[name],color=str(measurement_number/np.max(datasets)*0.7), linewidth=1)
        plt.plot(0,0,label='First',color=str(0.7/np.max(datasets)))
        plt.plot(0,0, label='Last',color=str(0.7))
        plt.axis([0,2*t_background+t_flag1_open,0,None])
        plt.legend(loc='upper right')
        plt.title(str(self.position))
        if save:
            if not os.path.exists(savefolder):
                os.makedirs(savefolder)
            plt.savefig(savefolder+'measurement_set_'+str(self.position)+'.png', dpi=500)
        plt.show()
        plt.close()
        
    def sum_measurements(self):
        self.sum = sum(list(self.datadict.values()))
        
    def flip_normalize(self):
        """      
         
        """        
        #calculate some indices
        open_flag1 = int(np.round(self.t_background / self.timestep))  
        num_base = int(np.round(open_flag1 * 0.8)) #number of datapoints for averaging base. Factor of 0.8 can be changed
        open_flag2 = int(np.round(self.t_open_flag2 / self.timestep))
        num_top = int(np.round(open_flag2 * 0.5)) #number of datapoints for averaging top. Factor of 0.7 can be changed
        num_shift = int(np.round(open_flag2 * 0.4))
        
        #subtract baseline
        left = np.average(self.sum[:num_base])
        right = np.average(self.sum[-num_base:])
        slope = (right-left)/(self.t_flag1_open/self.timestep)
        
        self.baseline = np.arange(len(self.sum))*slope - slope*open_flag1 + left        
        self.normalizeddata = self.sum-self.baseline
        
        #normalize
        norm = np.average(self.normalizeddata[open_flag1+num_shift:open_flag1+num_shift+num_top])
        self.normalizeddata /= norm
        
        #flip
        self.normalizeddata -= 1 #1 because normalized
        self.normalizeddata = -self.normalizeddata
        
        
    def fit_stickingprob(self, t_fit=1.8, t_shift=0.6):  
        num_fit = int(np.round(t_fit/self.timestep))
        num_shift = int(np.round(t_shift/self.timestep))
        index = int(np.round((self.t_background + self.t_open_flag2)/self.timestep)) + num_shift
        
        fitparams, cov= np.polyfit(list(self.timedict.values())[0][index:index+num_fit],self.normalizeddata[index:index+num_fit], 1, cov=True)
        slope = fitparams[0]
        y = fitparams[1]
        
        self.fit = list(self.timedict.values())[0]*slope+y               
        t = self.t_background+self.t_open_flag2
        
        self.stickingprob = t*slope+y      
        self.sigma = np.sqrt(t**2*cov[0,0]+cov[1,1]+t*2*cov[0,1]) 
        
        
    def plot_analyzed_data(self, save=False):
        plt.plot(list(self.timedict.values())[0], self.sum/np.max(self.sum),label='Original summed data (a.u.)')
        plt.plot(list(self.timedict.values())[0], self.normalizeddata, label='Processed data')
        plt.plot(list(self.timedict.values())[0], self.baseline/np.max(self.sum), label='Baseline')
        plt.plot(list(self.timedict.values())[0], self.fit, label = 'Sticking dip fit')
        plt.scatter(self.t_background+self.t_open_flag2, self.stickingprob, label = 'Initial sticking probability')
        plt.legend(loc='best')
        plt.axis([0,None,0,1])
        plt.title('Position = '+str(self.position)+ ', Sticking probability = '+str(np.round(self.stickingprob,3)))
        if save:
            if not os.path.exists(savefolder):
                os.makedirs(savefolder)
            plt.savefig(savefolder+'analyzed_'+str(self.position)+'.png', dpi=500)        
        plt.show()
        plt.close()
    
########## End classes ###################
    
########## Functions #####################
        
def read_data_files(measurements, numbers, filenamestart='KW'): 
    """
    t_background: number of seconds of background before and after each measurement
    returns: measurements, a dictionary with positions as keys, and measurement
    object as values
    """

    if numbers == 'odd':
        filenames = number_to_str(datasets[np.argwhere(datasets%2)],filenamestart=filenamestart) 
        positions = positions_odd
            
    if numbers == 'even':
        filenames = number_to_str(datasets[np.argwhere(datasets%2==0)],filenamestart=filenamestart) 
        positions = positions_even
     
    timestep = 0
    for filename in filenames:
        if not os.path.exists(folder+filename+ext):
            print(folder+filename+ext)
            continue
        time, data = np.loadtxt(folder+filename+ext, skiprows=3, unpack=True, usecols=(0,1)) #read file    
        index_start, index_stop = find_times_multiple(data) #find indices of rising and falling edge of the data
        print (filename)
        if not timestep: #ensures timestep is determined only once, with the first dataset.
            timedif = np.roll(time,-1)-time       
            timestep = np.min(timedif[:-1])


        for pos, index in zip(positions, index_start):
            measurements[pos].timestep = timestep
            times = np.arange(time[index]-measurements[pos].t_background, 
                              time[index]+measurements[pos].t_background+measurements[pos].t_flag1_open+0.000000001, 
                              timestep) #make time array for interpolation   
            measurements[pos].datadict[filename] = np.interp(times, time, data)
            measurements[pos].timedict[filename] = times - times[0]
            
        
def number_to_str(numbers, filenamestart=''):
    """
    numbers: array or list with numbers
    Converts list of numbers to list of filenames (in the standard form of KW01 etc)
    """
    strs = []
    for i in range(len(numbers)):
        name = filenamestart + str(int(numbers[i]/10)) + str(int(numbers[i]%10))
        strs.append(name)
    return strs

    
def find_times_multiple(data):
    """
    Output: an array of start and an array of end indices for each measurement.
    
    If necessary, the threshold for rising and falling edge can be determined independently (but now it isn't)
    """
    thresh = (np.max(data)+np.min(data))/4 #/4 because the higher part of the data is noisier than the lower part. and sticking dip    
    cond_up1 = data < thresh
    cond_up2 = np.roll(data,-1) > thresh
    index_up = np.argwhere(cond_up1 * cond_up2)#indices for rising edge. gives the index left of the transition
    
    cond_down1 = data > thresh
    cond_down2 = np.roll(data,-1) < thresh
    index_down = np.argwhere(cond_down1 * cond_down2) #index left of the transition, falling edge
    
    return np.array(index_up), np.array(index_down)

main()