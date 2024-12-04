# -*- coding: utf-8 -*-
"""
Created on Mon May 20 16:30:21 2019

@author: jansenc3
"""

import numpy as np
from matplotlib import pyplot as plt

folder = 'P:/DATA/2019/05 May/190516/KWTPD/'
folder = 'P:/DATA/2020/01 Jan/200113/KW/'
filename = 'KW01'
ext = '.txt'

t_open_flag1 = 1
t_open_flag2 = 3 #s, time after opening of flag 1
t_dose = 3 #s, dose time of the beam 



def main():
    time, data = np.loadtxt(folder+filename+ext, skiprows=3, unpack=True)
    
    open_flag1, close_flag1 = find_times_single(data)
    
    data_proc = flip_normalize(data, open_flag1, close_flag1)
    
    stickingprob, fit = fit_stickingprob(time, data_proc, open_flag1, close_flag1, t_open_flag2)
    
    index_start, index_end = find_times_multiple(data)
    

    
    plt.plot(time, data/np.max(data),label='Original data (a.u.)')
    plt.plot(time, data_proc, label='Processed data')
    plt.plot(time, fit, label = 'Sticking dip fit')
    plt.scatter(time[open_flag1]+t_open_flag2, stickingprob, label = 'Initial sticking probability')
    plt.scatter(time[index_start], index_start/(index_start*2), label = 'open flag 1')
    plt.scatter(time[index_end], index_end/(index_end*2), label = 'close flag 1')
    plt.legend(loc='best')
    plt.axis([None,60,0,1])
    plt.title('Sticking probability = '+str(np.round(stickingprob,3)))
    plt.show()
    plt.close()
    
def find_times_multiple(data):
    """
    Output: an array of start and a list of end indices for each measurement
    """
    thresh = (np.max(data)+np.min(data))/3 #/3 because the higher part of the data is noisier than the lower part    
    cond_up1 = data < thresh
    cond_up2 = np.roll(data,-1) > thresh
    index_up = np.argwhere(cond_up1 * cond_up2)#indices for rising edge. gives the index left of the transition
    
    cond_down1 = data > thresh
    cond_down2 = np.roll(data,-1) < thresh
    index_down = np.argwhere(cond_down1 * cond_down2) #index left of the transition, falling edge
    
    return np.array(index_up), np.array(index_down)

def find_times_single(data):
    """
    Only works for a single KW measurement. If there are multiple measurements, 
    they have to be separated first.
    threshold for detecting opening of flag: peak/2. 
    both indices are one datapoint to the left of the threshold, because 
    opening the flag happens before the change in partial pressure
    """
    thresh = (np.max(data)+np.min(data))/2
    ind = np.argwhere(data>thresh)
    open_flag1 = int(ind[0])-1 #-1 to get the index to the left of the transition
    close_flag1 = int(ind[-1])
    return open_flag1, close_flag1
    

def flip_normalize(data, open_flag1, close_flag1, num_avg=20, num_shift=5):
    """
    num_avg: number of points left and right of the pressure jump that
    are used for averaging the baseline level
    num_shift: determines how far left and right the averaging is done, to 
    prevent including part of the jump/drop
    """
    print('not finished')
    #subtract baseline
    left = np.average(data[open_flag1-num_shift-num_avg:open_flag1-num_shift])
    right = np.average(data[close_flag1+num_shift:close_flag1+num_shift+num_avg])
    
    slope = (right-left)/(close_flag1-open_flag1)
    
    baseline = np.arange(0,slope*len(data),slope) - slope*open_flag1 + left
    
    data = data-baseline
    
    #normalize
    norm = np.average(data[open_flag1+num_shift:open_flag1+num_shift+num_avg])
    data /= norm
    
    #flip
    data -= 1 #1 because normalized
    data = -data
    
    return data

def fit_stickingprob(time, data, open_flag1, close_flag1, t_open_flag2, num_fit=5, num_shift=5):
    index = np.argmax(data[open_flag1+num_shift:close_flag1-num_shift])
    index += open_flag1+num_shift
    
    slope, y = np.polyfit(time[index:index+num_fit],data[index:index+num_fit], 1)

    fit = time*slope+y
    
    stickingprob = y + slope*(time[open_flag1]+t_open_flag2)
    
    return stickingprob, fit
    
    
    
main()