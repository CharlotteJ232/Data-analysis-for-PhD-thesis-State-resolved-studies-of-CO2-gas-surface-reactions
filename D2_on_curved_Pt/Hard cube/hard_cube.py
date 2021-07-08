# -*- coding: utf-8 -*-
"""
Created on Fri May  1 14:39:17 2020

@author: Charlotte
"""
from sympy import *
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

#constants 
k = 1.38064852E-23
R = 8.31446261815324E-3 #kJ/K/mol
amu = 1.66053904E-27
avo = 6.02214E23

#parameters
ms_amu = 195 #amu
mm_amu = 4 #amu
Ts = 150 #Kelvin
#Em_kjm = 0.8 #kJ/mol
Ewell_kjm = 12.06 #kJ/mol
Em_list_kjm = np.array([0.7,0.8,1.2,2,2.9,4.0, 4.77,6.1])
#Em_list_kjm = np.array([0.8])


#Calculated values to SI units
ms = ms_amu * amu
mm = mm_amu * amu
#Em = Em_kjm / avo * 1000
Ewell = Ewell_kjm /avo * 1000
Em_list = Em_list_kjm / avo * 1000

Emax = 5*k*Ts #maximum energy in distribution



#speed distribution of surface
speed_range = np.linspace(-np.sqrt(2*Emax/ms),np.sqrt(2*Emax/ms),100)
#collision_prob = np.absolute((v0m+speed_range)/(2*v0m)) #because the molecule is more likely to hit the surface molecule while it is moving down
speed_dist = np.exp(-ms*speed_range**2/(2*k*Ts)) #speed distribution of the surface atoms
v1s_range = np.zeros(len(speed_range))
v1m_range = np.zeros(len(speed_range))


trapping_prob_list = []



def main():   
#    print ('speed range ',speed_range, 'speed dist ', speed_dist)

    for Em in Em_list:
        print (Em)
        v0m = -np.sqrt(2*(Em+Ewell)/mm)
        collision_prob = np.absolute((v0m+speed_range)/(2*v0m)) #because the molecule is more likely to hit the surface molecule while it is moving down

        for i in range(len(speed_range)):
            v0s = speed_range[i]
            v1m = Symbol('v1m')
            v1s = Symbol('v1s')
            pars = [v1m, v1s]
            eqs = [mm*v0m**2+ms*v0s**2-mm*v1m**2-ms*v1s**2,
                   mm*v0m+ms*v0s-mm*v1m-ms*v1s]
    
            
            solutions = solve(eqs, pars)
            
            v1m_range[i], v1s_range[i] = filter_trivial(solutions, v0m)
            
        trapped = 0.5*v1m_range**2*mm < Ewell
        
        trapping_prob = np.sum(trapped*collision_prob*speed_dist)/np.sum(collision_prob*speed_dist)

        
        trapping_prob_list.append(trapping_prob)
        
    print(trapping_prob_list)
    
    X = np.column_stack((Em_list_kjm, np.array(trapping_prob_list)))   
    np.savetxt(folderstart+folder+'trapping_prob_'+str(Ewell_kjm)+'_kjm.txt',X, header='E(kJ/mol), trapping probability')
    
    plt.scatter(Em_list_kjm, trapping_prob_list)
    plt.title('Well depth = '+str(Ewell_kjm)+' kJ/mol')
    plt.xlabel('Molecule kinetic energy (kJ/mol)')
    plt.ylabel('Trapping probability')
    plt.savefig(folderstart+folder+'trapping_prob_'+str(Ewell_kjm)+'_kjm.png')
    plt.show()
    plt.close()
    
    
def filter_trivial(solutions, v0m):
    """
    filters out the trivial solution where the final velocities are equal to the
    starting velocities
    """
    solutions_filtered=[]
    for sol in solutions:
        if np.absolute(sol[0]-v0m) < 1E-5:
            continue
        solutions_filtered.append(sol)
    if len(solutions_filtered)==1:
        return solutions_filtered[0]
    else:
        print('multiple solutions')
        return solutions_filtered

    
main()