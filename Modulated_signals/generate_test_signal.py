# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import numpy as np
from matplotlib import pyplot as plt
from scipy import signal as sn

measure_time = 5000 #s
measure_rate = 100 #samples/s
modulation_freq = 0.5 #hz
modulation_amplitude = 0.00001 #fraction of the total signal
noise_amplitude = 100 #times larger than the modulation amplitude

folderstart = 'D:/Surfdrive/DATA/'
savefolder = folderstart+'Testdata/modulated_kw/'

time = np.arange(0,measure_time,1/measure_rate)
signal = 1 - sn.square(time*modulation_freq*2*np.pi)*modulation_amplitude/2 - modulation_amplitude/2
noise = np.random.normal(scale=noise_amplitude*modulation_amplitude,size=len(time))

plt.plot(time, signal+noise)
plt.plot(time, signal)
plt.axis([0,20,0,np.max(signal+noise)*1.1])
plt.show()
plt.close()

savedata = np.column_stack((time, signal+noise, signal))
np.savetxt(savefolder+'rate_'+str(measure_rate)+'_freq_'+str(modulation_freq)+'_modulation_'+str(modulation_amplitude)+'_noise_'+str(noise_amplitude)+'.txt',
           savedata, header='time, noisy signal, clean signal')