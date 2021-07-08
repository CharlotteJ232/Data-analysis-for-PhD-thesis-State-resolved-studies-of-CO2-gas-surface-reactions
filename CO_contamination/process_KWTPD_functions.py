# -*- coding: utf-8 -*-
"""
Created on Mon Mar 25 15:32:09 2019

@author: jansenc3
"""

import numpy as np



def remove_overr(array):
    datamask = np.argwhere(array < 1e30) #to get all indices of valid data points
    datamask = datamask[:,0] #to get only row indices
    return array[datamask,:] #to return only rows where no nonsense data is present

def get_normalizing_factor(refdata):
    weight = 0.7
    refdata = refdata[refdata<1e30]#to remove any nonsense data
    maxpos = np.argmax(refdata)
    yvalue = weight*np.max(refdata)+(1-weight)*refdata[-1] #weighted average, values are a bit arbitrary
    pos = maxpos
    while refdata[pos] > yvalue:
        pos+=1
    
    array = np.zeros(3)
    for direction in [1,-1]:
        k=pos
        while True:
            left = refdata[k]-refdata[k-1]
            right = refdata[k+1]-refdata[k]
            if left*right < 0:
                array[direction]=refdata[k]
                break
            k += direction    
        
    factor = np.absolute(array[1]-array[-1])
    print ('factor ', factor)
    return factor

def convert_temperature(Tmin, Tmax, data):
    deltaT = Tmax - Tmin
    deltaV = np.max(data) - np.min(data)
    data = data/deltaV*deltaT
    data = data-np.min(data)+Tmin  
    return data
    
def define_peak(peak_pos,data,side='r'):
    """
    Very simple at this point.
    checks where left and right side of (local) dip are approximately equal.
    """
    
    if side == 'l':
        direction = -1
        factor = 1.5
    elif side == 'r':
        direction = 1
        factor = 1/1.5
    
    k = peak_pos 
    while True: #going right
        left = data[k-1]-data[k]
        right = data[k+1]-data[k]
        if left > 0:
            if right > 0:
                if direction*right > direction*factor*left: 
                    edge = k
                    break
        k += direction    

    return edge

def define_peak_smooth(peak):
    print ('not finished')
    
def find_peak_positions(data):
    print ('not finished')

def integral_linear_background(peak_pos, time, data, returnlr=False):
    """
    Draws a line between left and right edge of the peak, subtracts that line 
    from the data. then takes integral
    """   
    left = define_peak(peak_pos,data,side='l')
    right = define_peak(peak_pos,data,side='r')
    
    peak = data[left:right]
    peaktime = time[left:right]
    
    slope = (data[left]-data[right])/(left-right)
    line = np.arange(data[left],data[right],slope)
    
    peak = peak - line
    
    dtime = peaktime - np.roll(peaktime,1)
    dtime = dtime[1:] #because first element is bullshit
    
    dpeak = (peak + np.roll(peak,1))/2 #average 
    dpeak = dpeak[1:]
    
    integral = np.sum(dtime*dpeak)
    
    if returnlr:
        return integral, left, right
    else:    
        return integral
    


def remove_CO2_background(temp, CO2, data):
    #select part of the data to compare the two signals: the horizontal parts of the temperature signal
    #so left of the minimum and right of the maximum of the signal should be ok
    maxi = np.argmax(temp)
    mini = np.argmin(temp)
    #divide two signals at these parts
    compare_data = np.concatenate((CO2[:mini]/data[:mini],CO2[maxi:]/data[maxi:]))

    divide = np.average(compare_data)
    #divide CO2 signal by the average of this
    CO2 = CO2/divide
    #subtract from CO data
    data = data-CO2
    #correct for background coming from actual CO? 
    #subtract integral from signal? 
    return data
    