# -*- coding: utf-8 -*-
"""
Created on Mon May 20 16:30:21 2019

@author: jansenc3
"""

import numpy as np
from matplotlib import pyplot as plt
import os

folderstart = 'P:/Surfdrive/DATA/'
folderstart = 'C:/Users/Werk/surfdrive/DATA/'

folder = 'P:/DATA/2020/02 Feb/200207/KW/'
datasets = [9,10]

folder = folderstart+'2021/06 Jun/210618/KW/'
dip_depth = 4 #used in cutting the dataset into separate measurements
measurement_names = {'laseron':['KW09'], 'laseroff':['KW10']}
# measurement_names = {'laseron_indiv':['KW06','KW07','KW08'], 'laseroff_indiv':['KW10']}

folder = folderstart+'2021/06 Jun/210622/KW/'
dip_depth = 8
measurement_names = {'laseron':['KW02', 'KW03'], 'laseroff':['KW01','KW04']}

filenamestart = 'KW' #'KW' for standard 'KW01' filenames
ext = '.txt'
savefolder = folder+'Images/'

t_flag1_open = 9 #s, total time flag 1 is open
t_open_flag2 = 3 #s, time after opening of flag 1
#t_dose = 3 #s, dose time of the beam 
t_background = 4 #amount of seconds left and right of the measurement that are kept in the datasets. this can be chosen

save_individual = True



def main():
    if not os.path.exists(savefolder):
        os.makedirs(savefolder)    
    
    
    measurements = {}
    for measurement_name in measurement_names.keys():
        measurements[measurement_name] = measurement(measurement_name, t_background, t_flag1_open, t_open_flag2)
        plot_raw_data(measurement_names[measurement_name])
        read_data_files(measurements[measurement_name], measurement_names[measurement_name])

    stickinglist = []
    errorlist = []
    for measurem in measurements.values():
        measurem.plot_measurement_set(save=save_individual)    
        measurem.sum_measurements()
        measurem.flip_normalize()
        measurem.fit_stickingprob(t_fit=2.2,t_shift=0.5)
        measurem.plot_analyzed_data(save=save_individual)
        stickinglist.append(measurem.stickingprob)
        errorlist.append(measurem.sigma)

    for measurem in measurements.values():
        total_time = 2*t_background+t_flag1_open
        index = int((t_background+0.5*t_open_flag2)/total_time * len(measurem.sum))
        plotdata = measurem.sum-np.min(measurem.sum)
        plotdata /= np.average(plotdata[index-10:index+10])
        plt.plot(list(measurem.timedict.values())[0], plotdata, label=measurem.position)
    plt.axis([0,total_time,0,None])
    plt.xlabel('Time (s)')
    plt.ylabel('QMS signal (normalized)')
    plt.legend()
    plt.savefig(savefolder+'comparison.png',dpi=500)
    plt.show()
    plt.close()

    stickinglist = np.array(stickinglist) 
    
    plt.scatter([1,2], stickinglist)
    plt.errorbar([1,2], stickinglist, yerr=errorlist, fmt='none', capsize=5)
    plt.axis([None, None, 0, None])
    plt.title('sticking probability as a function of crystal position')
    plt.savefig(savefolder+'stickingprob_vs_position.png', dpi=500)
    
    # sav = np.column_stack((positions, stickinglist, errorlist))
    # np.savetxt(savefolder+'stickingprob_vs_position.txt',sav, header='Position (mm or mm/10), sticking probability, error')
    plt.show()
    plt.close()   

########## End main ######################
    
########## Classes  ######################
    
class measurement:
    """
    Contains all the information for a measurement with certain experimental conditions
    Can contain several datasets, which can then be added together for averaging
    """
    def __init__(self, position, t_background, t_flag1_open, t_open_flag2):
        self.position = position
        self.t_background = t_background
        self.t_flag1_open = t_flag1_open
        self.t_open_flag2 = t_open_flag2

        self.timedict = {}
        self.datadict = {}
  
    
    def plot_measurement_set(self, save=False):
        col = 0 #grayscale color
        for name in self.datadict.keys():
            measurement_number = float(name.replace(filenamestart,''))
            left = self.timedict[name]<t_background*0.9 #*0.9 just to be sure I don't include a part of the measurement   
            baseleft = np.average(self.datadict[name][left]) 
            plt.plot(self.timedict[name],self.datadict[name]-baseleft, linewidth=0.5, color=str(col/len(self.datadict.keys())))
            col+=1
        plt.plot(0,0,label='First',color=str(0/len(self.datadict.keys())))
        plt.plot(0,0, label='Last',color=str((col-1)/len(self.datadict.keys())))
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

def plot_raw_data(filenames):
    for filename in filenames: 
        time, data = np.loadtxt(folder+filename+ext, skiprows=3, unpack=True, usecols=(0,1)) #read file    
        plt.plot(time,data)
    plt.title(filenames)
    plt.show()
    plt.close()
        
def read_data_files(measurementname, filenames): 
    """
    t_background: number of seconds of background before and after each measurement
    returns: measurements, a dictionary with positions as keys, and measurement
    object as values
    """
     
    timestep = 0
    for filename in filenames:
        time, data = np.loadtxt(folder+filename+ext, skiprows=3, unpack=True, usecols=(0,1)) #read file    
        index_start, index_stop = find_times_multiple(time, data) #find indices of rising and falling edge of the data
        if not timestep: #ensures timestep is determined only once, with the first dataset.
            timedif = np.roll(time,-1)-time       
            timestep = np.min(timedif[:-1])


        for i in range(len(index_start)):
            measurementname.timestep = timestep
            times = np.arange(time[index_start[i]]-measurementname.t_background, 
                              time[index_start[i]]+measurementname.t_background+measurementname.t_flag1_open+0.000000001, 
                              timestep) #make time array for interpolation   
            measurementname.datadict[filename+str(i)] = np.interp(times, time, data)
            measurementname.timedict[filename+str(i)] = times - times[0]
            

        
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

    
def find_times_multiple(time, rawdata):
    """
    threshold = np.max(data)/threshvalue
    Output: an array of start and an array of end indices for each measurement.
    
    If necessary, the threshold for rising and falling edge can be determined independently (but now it isn't)
    """
    baseline = rawdata[0] + (rawdata[-1]-rawdata[0])/(time[-1]-time[0]) * time
    data = rawdata - baseline
    thresh = np.max(data)/dip_depth #min counts stronger because the higher part of the data is noisier than the lower part. and sticking dip    
    
    cond_up1 = data < thresh
    cond_up2 = np.roll(data,-1) > thresh
    index_up = np.argwhere(cond_up1 * cond_up2)#indices for rising edge. gives the index left of the transition
    
    cond_down1 = data > thresh
    cond_down2 = np.roll(data,-1) < thresh
    index_down = np.argwhere(cond_down1 * cond_down2) #index left of the transition, falling edge
    
    return np.array(index_up), np.array(index_down)

main()