# -*- coding: utf-8 -*-
"""
Created on Fri Apr 24 11:25:08 2020

@author: Charlotte
"""

import numpy as np

E_beam = 80 #meV
well_depth = 65 #meV
E_in = E_beam+well_depth

def main():
    vm, vs, theta_m, theta_s = np.loadtxt('solutionslist.txt',unpack=True, skiprows=1)
    dtheta = theta_s[1]-theta_s[0]
    
    area = np.absolute(np.sin(dtheta) * 2 * np.pi * np.sin(theta_s))
    print(area)
    
    E_out = (np.cos(theta_m)*vm)**2 * E_in
    
    print(E_out)
    
    trapped = E_out < well_depth
    
    print (trapped)
    
    trapping_prob = np.sum(trapped * area) / np.sum(area)
    
    print (trapping_prob)
    
    
main()