# -*- coding: utf-8 -*-
"""
Created on Wed Apr 24 13:50:45 2019

@author: jansenc3
"""
import glob
import numpy as np
from matplotlib import pyplot as plt


files = ['01','02','03','04','05','06']
positions = [46, 35, 25, 15, 5, 0]
number_of_cycles = 1.875 #number of chopper cycles in dataset. Usually 1.875.
number_of_peaks = number_of_cycles*4

smallpeakwidth = 0.2 #peak bottom width in ms. use None if you don't know.
largepeakwidth = None


######### Main script ########################################################

def main():
    time, data = read_files('P:/DATA/2019/04 Apr/190423',files=files)
    peak_index = find_peaks(data, number_of_peaks)  
    
    
    plt.plot(time[0], data[0], c='r')
    plt.scatter(time[0][peak_index[0]],data[0][peak_index[0]])
    plt.axis([0.2,0.5,500,1000])
    plt.show()
    plt.close()
   
########## Functions #########################################################

def read_files(folder, filename='TOF', files=None, ext='.asc'):
    paths = []
    time = [] #This will be a list of arrays
    data = [] #List of arrays
    
    if files==None:
        paths = glob.glob(folder+'/'+filename+'*'+ext)
    else:
        for i in range(len(files)):
            paths.append(folder+'/'+filename+files[i]+ext)
    
    for k in range(len(paths)):
        timetemp, datatemp = np.loadtxt(paths[k],unpack=True)
        time.append(timetemp)
        data.append(datatemp)
    
    return time, data

def find_peaks(data, number_of_peaks):
    peak_index = [] #list of arrays. each array contains the indices of all the peaks in one dataset
    
    for dataset in data:
        window = int(len(dataset)/number_of_peaks) #calculate window size for a single peak
        peaks = np.empty(int(np.ceil(number_of_peaks)))
        for i in range(int(np.ceil(number_of_peaks))):
            peaks[i] = window*i + np.argmax(dataset[window*i:window*(i+1)]) #searches peak in the right window and saves index
        peak_index.append(np.array(peaks.astype(int)))
    
    return peak_index

def fit_peaks(data, peak_index, number_of_peaks):
    print ('not finished')
    
def peakshape():
    print ('not finished')

if __name__ == "__main__":
    main()