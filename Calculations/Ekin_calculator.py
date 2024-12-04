import numpy as np
from matplotlib import pyplot as plt

R = 8.31446261815324

#factor is flow factor of MFC, degf is number of degrees of freedom that can be cooled
gas_info = {'N2':   {'mass':28, 'factor':1,     'degf':7},
            'CO2':  {'mass':44, 'factor':0.74,  'degf':7},
            'He':   {'mass':4,  'factor':1.386, 'degf':5},
            'O2':   {'mass':32, 'factor':0.988, 'degf':7},
            'H2':   {'mass':2,  'factor':1.008, 'degf':7},
            'D2':   {'mass':4,  'factor':0.995, 'degf':7},
            'Ar':   {'mass':40, 'factor':1.395, 'degf':5},
            'CO':   {'mass':28, 'factor':0.995, 'degf':7}}

##### Parameters #####

Tnozzle = 1050
Tnozzle = 850
Tnozzle = 600
energy_efficiency = 0.9 #supplementary fig 3-5 of nakamura paper

#For each gas in the beam, put in an array the different flows 
#you want to calculate. The code will calculate the energy for
#the fist set of flows, then for the second set, etc. 
beam = {'CO2':  np.array([7, 3, 2, 1, 1, 1, 1]),
        'He':   np.array([12,15, 15,30, 8, 12, 14])}

beam = {'CO2':  np.array([0.5, 0.5, 0.5, 0.5, 0.7, 0.5, 0.3]),
        'He':   np.array([7, 8.64, 26.5, 6.4, 7.1, 9, 12])}

# beam = {'O2': np.array([15, 10, 7.5, 5, 2, 2]),
#         'Ar': np.array([0, 5, 7.5, 10, 20, 5])}

beam = {'D2': np.array([10]),
        'D2': np.array([10])}

beam = {'CO2':  np.array([5]),
        'He':   np.array([9])}
 
gas1 = 'O2'
gas2 = 'Ar'

gas1 = 'CO2'
gas2 = 'He'

# gas1 = 'D2'
# gas2 = 'D2'

##### Calculations #####

m1 = gas_info[gas1]['mass']
flow1 = beam[gas1]*gas_info[gas1]['factor']
degf1 = gas_info[gas1]['degf']

m2 = gas_info[gas2]['mass']
flow2 = beam[gas2]*gas_info[gas2]['factor']
degf2 = gas_info[gas2]['degf']

#### calculations ####

x1 = flow1/(flow1+flow2)
x2 = 1-x1
degf = degf1*x1 + degf2*x2

E = (m1/(x1*m1+x2*m2))*degf/2*R*Tnozzle
E_ev = E*0.01036410/1000
E_ev_adj = E_ev * energy_efficiency
velocity = np.sqrt(2*E*1000*energy_efficiency/m1) #in m/s

for i in range(len(flow1)):
        print(gas1,beam[gas1][i],'ml/min,', 
              gas2,beam[gas2][i],'ml/min,', 
              'E=',np.round(E[i]/1000,2), 'kJ/mol or',np.round(E_ev[i],2), 'eV,',
              'Adjusting for energy conversion efficiency E=', str(np.round(E_ev_adj[i],2)), 'eV',
              str(np.round(velocity[i],0))+' m/s',
              str(np.round(x1[i]*100,2))+'%',gas1,'in',gas2)

plt.plot(x1*100, E_ev_adj, 'o')
plt.xlabel('Percentage of '+gas1+' in '+gas2)
plt.ylabel('Energy adjusted for conversion efficiency (eV)')
plt.title(str(Tnozzle)+' K')
plt.show()
plt.close()


