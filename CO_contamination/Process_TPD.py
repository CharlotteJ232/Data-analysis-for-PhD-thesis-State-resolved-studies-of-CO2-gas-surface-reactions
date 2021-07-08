# -*- coding: utf-8 -*-
"""
Created on Tue Mar 19 11:34:14 2019

@author: jansenc3
"""

import numpy as np
from matplotlib import pyplot as plt
from scipy.signal import savgol_filter

import plots as ps
import process_KWTPD_functions as pr

folder = '0325/TPD/'
filenames = ['06','09',11,12,13,14]
filenames = ['02','03','04']
integralx = [1,5,10,20,40,0.17] #beam time in minutes for the files
integralx = [30,60,180]
Tmin = 200
Tmax = 704
column = 2
plot_all = True

integrals = np.zeros(len(integralx))


###### Start of script #######
for k,j in zip(filenames,range(len(filenames))):
    filename = str(k)+'.txt'
    print (filename)
    start_time = np.loadtxt(folder+'tpd'+filename, skiprows=1, max_rows=1)
    names = np.loadtxt(folder+'tpd'+filename,skiprows=2, max_rows=1, dtype=str, delimiter=',')
    data = np.loadtxt(folder+'tpd'+filename, skiprows=3)
    
    #data = pr.remove_overr(data) #this causes errors
  
#    #Scale data using KW data
#    kwdata = np.loadtxt(folder+'kw'+filename, skiprows=3)
#    factor = pr.get_normalizing_factor(kwdata[:,1])
#    data[:,2:] = data[:,2:]/factor
#    
#    #Remove CO2 background
#    data[:,3] = pr.remove_CO2_background(data[:,1],data[:,2],data[:,3])

    #Calculate integral of peak
    integrals[j], left, right = pr.integral_linear_background(
                                np.argmax(data[:,column]), 
                                data[:,0],data[:,column],returnlr=True)
    print ('Integral ',integrals[j])

    #Convert voltage to temperature for plotting
    temp = pr.convert_temperature(Tmin, Tmax, data[:,1])
    sort = np.argsort(temp)
    
### PLOT ###  
    plt.plot(temp, data[:,column],label=str(k))

#    #Plot vs time with peak edge points 
#    if plot_all:  
#        ps.plot_peak_edges(data[:,0], data[:,column], left, right, names[column], 
#                        title='TPD'+str(k)+' Integral = ' + str(np.round(integrals[j],2)), 
#                        save=True, filepath=folder+'Images/tpd'+str(k)+'_'+names[column]+'.png')
#   
#    #Plot vs temperature
#    ps.plot_vs_temp(temp[sort], data[sort,column],names[column],'TPD'+str(k),
#                save=True,filepath=folder+'Images/tpd'+str(k)+'_'+names[column]+'_vsTemp.png') 
#    
#    #test smooth data
#    data[:,column] = savgol_filter(data[:,column], 35, 4)
#    
#    ps.plot_vs_temp(temp[sort], data[sort,column],names[column],'TPD'+str(k),
#                save=False,filepath=folder+'Images/tpd'+str(k)+'_vsTemp.png')
#    
#    ps.plot_vs_temp(temp[sort], (data[sort,column]-np.roll(data[sort,column],-1)),names[column],'TPD'+str(k),
#                save=False,filepath=folder+'Images/tpd'+str(k)+'_vsTemp.png')
plt.legend()
plt.show()
plt.close()
ps.plot_integral_vs_beamtime(integralx,integrals,filepath=folder+'Images/Integrals_'+names[column]+'_vsTemp.png')    
