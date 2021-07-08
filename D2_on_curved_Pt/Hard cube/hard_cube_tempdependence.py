# -*- coding: utf-8 -*-
"""
Created on Fri May  1 14:39:17 2020

@author: Charlotte
"""

import numpy as np
from matplotlib import pyplot as plt

"""
NOTES
Well depth terrace 30-40 meV2 - 3.38 kJ/mol
120-130 bottom step - 12.06 kJ/mol
170-180 top step - 16.88 kJ/mol
[3.38, 12.06, 16.88]

"""

#Data folders
folderstart = 'P:/DATA/'
folderstart = 'D:/Surfdrive/DATA/'
folder = 'D2 sticking probability as a function of step density and kinetic energy/Models/Hard cube/'
year = '2020/'

datadic = {0.7:{'kw_day':'03 Mar/200310', 'TOF_day':'03 Mar/200310','c':'red'},
           0.8:{'kw_day':'05 May/200512', 'TOF_day':'01 Jan/200128','c':'orange'},
           1.2:{'kw_day':'03 Mar/200305', 'TOF_day':'02 Feb/200220','c':'yellow'},
           2.0:{'kw_day':'06 Jun/200602', 'TOF_day':'05 May/200511','c':'green'},
           2.9:{'kw_day':'03 Mar/200303', 'TOF_day':'02 Feb/200218','c':'blue'},
           4.0:{'kw_day':'03 Mar/200311', 'TOF_day':'03 Mar/200306','c':'indigo'},
           4.77:{'kw_day':'05 May/200519', 'TOF_day':'05 May/200519','c':'purple'},
           6.1:{'kw_day':'05 May/200520', 'TOF_day':'02 Feb/200204','c':'violet'}}

#constants 
k = 1.38064852E-23
R = 8.31446261815324E-3 #kJ/K/mol
amu = 1.66053904E-27 
avo = 6.02214E23

#parameters
ms_amu = 195 #amu
mm_amu = 4 #amu
# Ts = 150 #Kelvin
Ewell_kjm = 6 #kJ/mol
Em_list_kjm = np.array([0.7,0.8,1.2,2,2.9,4.0, 4.77,6.1])


#Calculated values to SI units
ms = ms_amu * amu
mm = mm_amu * amu
Ewell = Ewell_kjm /avo * 1000
Em_list = Em_list_kjm / avo * 1000



trapping_prob_list = []

averages=False #set to true to show how averaged surface and molecule speeds do not work
Eavg=0.8
Ts_list=np.linspace(100,350,6)

def main():       
    for Ts in Ts_list:
        print (Ts)
        Emax = 5*k*Ts #maximum energy in surface distribution. the factor 5 is a choice I made
        #speed distribution of surface
        v0s_range = np.linspace(-np.sqrt(2*Emax/ms),np.sqrt(2*Emax/ms),100) #all surface atom speeds in an array. includes positive and negative
        # v0s_range = np.array(-np.sqrt(2*Emax/ms),np.sqrt(2*Emax/ms))
        v0s_dist0 = np.exp(-ms*v0s_range**2/(2*k*Ts)) #speed distribution of the surface atoms. The total does not equal 1

        
        Em_kjm, Em_dist = np.loadtxt(folderstart+year+datadic[Eavg]['TOF_day']+'/TOF/Images/4.0/energy.txt',unpack=True, skiprows=1)
        Em = Em_kjm/avo*1000
        
        v0m_range = -np.sqrt(2*(Em+Ewell)/mm) #initial velocity of molecule, defined as negative 
       
        v0m, v0s = np.meshgrid(v0m_range, v0s_range) #make 2D arrays of all combinations of surface and gas speeds
        v0m_dist, v0s_dist = np.meshgrid(Em_dist, v0s_dist0) #2D arrays of distributions        
        
        if averages:
            v0m = np.average(v0m)
            v0m_dist = np.average(v0m_dist)
            v0s = np.average(v0s)
            v0s_dist = np.average(v0s_dist)
        
        collision_prob = np.absolute((v0m-v0s)/(2*v0m)) #because the molecule is more likely to hit the surface molecule while they are moving towards each other        
        print(collision_prob)
       
        v1m = (mm - ms)/(mm+ms)*v0m + 2*ms/(mm+ms)*v0s #formula from wikipedia               
        trapped = 0.5*v1m**2*mm < Ewell #trapped if final energy is smaller than well depth
        
        trapping_prob = np.sum(trapped*collision_prob*v0m_dist*v0s_dist)/np.sum(collision_prob*v0m_dist*v0s_dist) #include energy distributions and probability 
        trapping_prob_list.append(trapping_prob)
        
    # print(trapping_prob_list)

    X = np.column_stack((Ts_list, np.array(trapping_prob_list)))   
    np.savetxt(folderstart+folder+'trapping_prob_'+str(Ewell_kjm)+'_kjm_temperature.txt',X, header='E(kJ/mol), trapping probability')
    
    
    #plot
    plt.scatter(Ts_list, trapping_prob_list)
    plt.title('Well depth = '+str(Ewell_kjm)+' kJ/mol')
    plt.xlabel('Surface temperature')
    plt.ylabel('Trapping probability')
    plt.savefig(folderstart+folder+'trapping_prob_'+str(Ewell_kjm)+'_kjm_'+str(Eavg)+'kjm_temperature.png')
    plt.show()
    plt.close()
    
    

    
main()