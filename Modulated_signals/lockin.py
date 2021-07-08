"""
Issues:
- Lock-in returns very high values when there is a relatively large dc component in the data. This should not happen. 
  May be a python issue

"""


import numpy as np
from matplotlib import pyplot as plt
import scipy.signal as sn
import os

folderstart = 'P:/Surfdrive/DATA/'
#folderstart = 'D:/Surfdrive/DATA/'
loadfolder = folderstart+'Testdata/modulated_kw/'
show_individual_plots = True
generate_override = False

measure_time = 50000 #s
measure_rate = 8 #samples/s
modulation_freq = 0.5 #hz
modulation_amplitude = 0.00001 #fraction of the total signal
noise_amplitude = 10 #times larger than the modulation amplitude

filepath = loadfolder+'rate_'+str(measure_rate)+'_freq_'+str(modulation_freq)+'_modulation_'+ str(modulation_amplitude)+'_noise_'+str(noise_amplitude)+'.txt'

def main():
  #generate data, if it does not exit yet
  if generate_override or not os.path.exists(filepath):
    print('Creating data..')
    time = np.arange(0,measure_time,1/measure_rate)
    signal = 1 - sn.square(time*modulation_freq*2*np.pi)*modulation_amplitude/2 - modulation_amplitude/2
    noise = np.random.normal(scale=noise_amplitude*modulation_amplitude,size=len(time))
    savedata = np.column_stack((time, signal+noise, signal))
    np.savetxt(filepath, savedata, header='time, noisy signal, clean signal')

  #load data
  time, noisy, clean = np.loadtxt(filepath,skiprows=1,unpack=True)
  noisy -= np.average(noisy) #because python has a problem with small relative variations
  clean -= np.average(clean)
  sine_for_testing = modulation_amplitude *  np.sin(time*modulation_freq*2*np.pi+np.pi) 

  number_of_windows = 10 #calculates everything 10 times, then averages
  start_times = np.arange(0, len(time), len(time)/number_of_windows).astype(int)
  a = np.flip(np.arange(10))
  integration_times = (len(time)/number_of_windows/1.5**a-1).astype(int) #1.5 is just a number that generated a nice distribution of integration_times
  correctness = np.zeros(len(integration_times))
  correctness_error = np.zeros(len(integration_times))
  correctness_clean = np.zeros(len(integration_times))
  correctness_sin = np.zeros(len(integration_times)) 

  #loop over different integration_times
  for j in range(len(integration_times)): 
    values_noisy = np.zeros(number_of_windows)
    for i in range(len(start_times)):
      #do lock-in
      values_noisy[i] = lockin(time, noisy, integration_times[j], start_time=start_times[i])
    #get average and std from all values
    value_noisy = np.average(values_noisy)
    error_noisy = np.std(values_noisy)
    
    #do lock-in
    value_clean = lockin(time, clean, integration_times[j])
    value_sin = lockin(time, sine_for_testing, integration_times[j])
    # print(value_noisy, value_clean, value_sin)

    # """
    # reference sine has an amplitude of 1
    # signal has an amplitude of modulation_amplitude/2
    # because it is a square wave, the sine in the signal has an amplitude of modulation_amplitude/2*4/pi
    # the expected output for 2 sine waves is Uout=1/2*Vsig*cos(phasediff)
    # so if we want Uout to be equal to modulation_amplitude, we need a conversion factor
    # Vsig = modulation_amplitude/2
    # Uout = modulation_amplitude = 1/2 * modulation_amplitude/2*4/pi * conversion_factor
    # conversion_factor = pi
    # """
    conversion_factor = np.pi

    value_noisy *= conversion_factor
    error_noisy *= conversion_factor
    value_clean *= conversion_factor
    value_sin *= 2

    correctness[j] = value_noisy/modulation_amplitude
    correctness_error[j] = error_noisy/modulation_amplitude
    correctness_clean[j] = value_clean/modulation_amplitude
    correctness_sin[j] = value_sin/modulation_amplitude


  #plot correctness as a function of integration time
  plt.errorbar(time[integration_times]/60, correctness, yerr=correctness_error, ls='o', capsize=5)
  plt.plot(time[integration_times]/60, correctness_clean, 'o')
  plt.plot(time[integration_times]/60, correctness_sin, 'o')
  plt.plot([0,np.max(time/60/number_of_windows)],[1,1],'--')
  plt.axis([0,np.max(time/60/number_of_windows),0,max(1.2,min(100,np.max(correctness)))])
  plt.title('S/N='+str(1/noise_amplitude)+', Modulation frequency='+str(modulation_freq)+' Hz, Sampling rate='+str(measure_rate)+' Hz')
  plt.xlabel('Integration time (minutes)')
  plt.ylabel('Lock-in result/modulation amplitude')
  plt.show()
  plt.close()



def lockin(time, dataset, integration_time, start_time=0):
    end_time = start_time + integration_time
    #generate sin and cos with correct frequency
    sinewave = np.sin(time*modulation_freq*2*np.pi)
    coswave = np.cos(time*modulation_freq*2*np.pi)

    val_sin = np.sum(dataset[start_time:end_time]*sinewave[start_time:end_time]) / len(time[start_time:end_time])
    val_cos = np.sum(dataset[start_time:end_time]*coswave[start_time:end_time]) / len(time[start_time:end_time])

    return np.sqrt(val_sin**2 + val_cos**2)
  

def plot_data(time, noisy, clean, sinewave):
  plt.plot(time, noisy+0.5)
  plt.plot(time, clean)
  plt.plot(time, 0.5-modulation_amplitude*0.5+modulation_amplitude*sinewave)
  # plt.axis([0,100,0,np.max(noisy)*1.1])
  plt.axis([0,20,None,None])
  plt.show()
  plt.close()

main()