# -*- coding: utf-8 -*-
"""
Created on Tue Feb  9 15:11:42 2021

@author: jansenc3
"""

import numpy as np

slaagkans = [0.48, 0.63]
faal_aantal = np.zeros(2)

for i in range(2):
    faalkans = 1-slaagkans[i]
    for j in range(51):
        kans = slaagkans[i]*faalkans**j
        faal_aantal[i] += j*kans 

print('Charlotte, Freek')
print('slaagkans ',slaagkans)
print('Faal_aantal ',faal_aantal)
print('Verhouding neg traits ',faal_aantal[0]/faal_aantal[1])
