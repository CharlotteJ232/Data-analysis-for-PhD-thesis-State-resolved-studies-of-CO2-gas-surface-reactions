# -*- coding: utf-8 -*-
"""

"""

import numpy as np
from matplotlib import pyplot as plt

approxfreq = 250

folders = ['P:/DATA/2019/05 May/190502/TOF/',
           'P:/DATA/2019/05 May/190517/Chopper/']
filenames = ['total', '3.32.55 PM']
labels = ['Closed box','Open box and fan']
ext = '.txt'

def main():
#    actualfrequency()
    plt.plot()
    plot()



def actualfrequency(): 
    for folder, filename in zip(folders, filenames):
        time, freq = np.loadtxt(folder+filename+ext, unpack=True)
        
        threshold = approxfreq*3
        
        freq[freq>threshold] /= 2
        freq /= 2
        
        data = np.zeros((len(time),2))
        data[:,0] = time
        data[:,1] = freq
        
        np.savetxt(folder+filename+'_actualfrequency'+ext, data)

#relfreq = freq / np.min(freq)

def plot():
    for folder, filename, label in zip(folders, filenames, labels):
        time, freq = np.loadtxt(folder+filename+'_actualfrequency'+ext, unpack=True)
        time -= time[0]
        time /= 3600
#        freq /= freq[int(len(freq)/2)]
        slope, y0 = fitline(time, freq)
        plt.plot(time, slope*time+y0, color='black')
        plt.plot(time, freq, label=label)
    plt.legend()
    plt.title('Chopper frequency')
    plt.axis([None, None, 254, 255])
    plt.xlabel('Time (h)')
    plt.ylabel('Relative frequency')
    plt.savefig(folder+'Images/chopperfreq.png',dpi=1000)
    plt.show()
    plt.close()
    
def fitline(time, freq):
    time = time[int(len(time)/2):]
#    time -= time[0]
    freq = freq[int(len(freq)/2):]
#        freq /= freq[0]
    slope, y0 = np.polyfit(time, freq, 1)
    
    print ('shift', slope/freq[0]*100, '% per minute')
    
    return slope, y0

main()






