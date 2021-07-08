# -*- coding: utf-8 -*-
"""
Created on Fri Jul  3 14:57:09 2020

@author: Charlotte
"""
import numpy as np
from matplotlib import pyplot as plt
plt.style.reload_library()
plt.style.use('voorbeeld')

eta = 1/1.35 #some dimensionless parameter, used in calculating tau
E_mev = 73 #in meV, used in calculating tau
Ld_v = 1E-14 #kortere afstanden -> effect van lagere s0 bij hogere T is minder sterk
step_density = np.linspace(0.00001, 0.5, 50) #nm-1
Ld = 1/step_density * 1E-9 #m
v = 1781 #m/s

#constants 
k = 1.38064852E-23
R = 8.31446261815324E-3 #kJ/K/mol
amu = 1.66053904E-27 
avo = 6.02214E23
h = 6.62607004E-34
mev_to_kjm = 96.487/1000

E = mev_to_kjm * E_mev/avo * 1000

Ts_list, S_trap = np.loadtxt('trapping_prob.txt', skiprows=1,unpack=True)

# plt.plot(Ts_list, S_trap, 'o', label='Trapping probability')

# tau = np.exp(eta*E/k/Ts_list)/k/Ts_list*h

# print(tau)

# s0 = (1-np.exp(-Ld_v/tau)) / Ld_v * tau

# print(s0)

# # plt.plot(Ts_list, s0)
# plt.plot(Ts_list, tau/np.max(tau), label='tau')
# plt.plot(Ts_list, s0, label='s0')

# plt.legend()
# plt.show()
# plt.close()


for Ts, S in zip(Ts_list, S_trap):
    tau = np.exp(eta*E/k/Ts)/k/Ts*h
    s0 = (1-np.exp(-Ld/(v*tau))) / Ld *( v * tau)
    plt.plot(step_density, s0, label='Ts = '+str(Ts))
plt.ylabel('probability of finding defect')
plt.xlabel('step density')
plt.legend()
    

