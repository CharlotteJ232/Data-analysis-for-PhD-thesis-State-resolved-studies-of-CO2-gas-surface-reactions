# -*- coding: utf-8 -*-
"""
Created on Thu May  9 14:23:47 2019

@author: Charlotte
"""
import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit

plotpeakshape = True

def main():
########### Settings ########################################################     
    
    data_folder = 'P:/DATA/2019/04 Apr/190430/TOF/'
    data_folder = 'P:/DATA/2019/07 Jul/190725/TOF/'
    data_folder = 'P:/DATA/2019/08 Aug/190808/TOF/'
#    data_folder = 'P:/DATA/2019/08 Aug/190812/TOF/'
#    data_folder = 'P:/DATA/2019/05 May/190502/TOF/'
#    data_folder = 'C:/Users/Charlotte/surfdrive/DATA/TOF testdata/'   
    nr_of_datasets = 14
    start_dataset = 1
    

    data_nrs = np.linspace(start_dataset, nr_of_datasets+start_dataset-1, num=nr_of_datasets) #numbers of datasets
    data_positions = np.full(nr_of_datasets,46)   
#    data_positions = [46, 35, 25, 15, 5, 0] #positions on moving stage (mm), corresponding to number of dataset above
    
    measured_chopper_freq = False
    len_dataset = 0.375 #number of chopper cycles. usually x.875
    chopper_freq = 242 #Hz
    gas_mass = 0.04 #kg/mol
    avg_mass = (4*0.004+2*0.04)/6
#    avg_mass = 0.004
    nozzle_temp = 298 #K
    t0_guess = 0.84 #one of the fit parameters that cannot be estimated easily 
    
#    data_vext = np.linspace(31.6, 9.6, num=nr_of_datasets)
#    data_vext = [7.79, 9.53, 11.78, 13.71, 15.69, 17.52, 19.66, 21.75, 23.74, 25.84, 27.50, 29.67, 31.58]
#    data_vext = [11.78, 13.71, 15.69, 17.52, 19.66, 21.75, 23.74]     
#    data_vext = [31.77, 29.39, 27.47, 25.54, 23.73, 21.80, 19.48, 17.40, 15.79, 13.63, 11.58, 9.63, 7.82]   #25 jul          
    data_vext = [7.8, 8.35, 9.45, 10.32, 11.61, 13.47, 15.27, 17.67, 19.90, 21.53, 23.44, 25.76, 27.85, 29.50] #8aug
#    data_vext = [31.84, 41.5, 52.3, 61.4, 70.5, 81.0, 1.338] #8 aug, vfoc
#    data_vext = [1.479, 14.2, 28.92, 43.5, 61.7, 79.3] #8 aug, vfoc
#    data_vext = [7.85, 8.47, 9.71, 10.57, 11.65, 13.73, 15.59, 17.90, 19.71, 21.76, 23.30, 25.61, 27.15, 30.16, 32.00] #12 aug, vfoc = 79.4
#    data_vext = [7.82, 8.34, 9.46, 10.67, 11.32, 13.35, 15.60, 17.24, 19.64, 21.53, 23.53, 25.43, 27.61, 29.00, 32.14] #12 aug, vfoc = 0.002
#    data_vext = [7.95, 8.52, 9.74, 10.75, 11.49, 13.32, 15.96, 32.19] #12 aug, velec = 32
#    data_vext = np.full(nr_of_datasets, 15)
    data_vext = np.array(data_vext) + 18.78


############# Script ########################################################

    testplot_fitfunctions()

    
#    dataid = Dataidentifier(data_nrs, data_positions, data_vext)
#    
#    exp = Experiment(data_folder, dataid, len_dataset, chopper_freq, gas_mass, avg_mass, nozzle_temp, t0_guess, measured_chopper_freq)
#    exp.plot_raw_data()

#    
#    for nr, dataset in exp.datadict.items():
#        print (nr)
#        dataset.makepeaks(exp)
#        dataset.fitpeaks(exp)
#    plt.legend()
#    plt.axis([0.6,0.64, 0.95,1])
#    plt.show()
#    plt.close()
     
##    exp.calc_v()   
#    exp.calc_Tion()
#    testplot_tion(np.arange(7, 32, 2))
#    exp.calc_Tslit()
    
#    testplot_chopper(exp)
    
################# Classes ####################################################


class Dataidentifier:
    def __init__(self, data_nrs, data_positions, data_vext, filenamestart='TOF', chopperfilename = 'Chopper', ext='.asc'):
        self.data_nrs = data_nrs
        self.data_positions = data_positions
        self.data_vext = data_vext
        self.filenamestart = filenamestart
        self.chopperfilename = chopperfilename
        self.ext = ext



class Experiment:
    def __init__(self, data_folder, dataid, len_dataset, chopper_freq, 
                 gas_mass, avg_mass, nozzle_temp, t0_guess, measured_chopper_freq):
        self.len_dataset = len_dataset
        self.chopper_freq = chopper_freq
        self.gas_mass = gas_mass
        self.avg_mass = avg_mass
        self.nozzle_temp = nozzle_temp
        self.t0_guess = t0_guess
        self.data_folder = data_folder
        self.measured_chopper_freq = measured_chopper_freq
        
        self.read_datasets(dataid)
        
    def read_datasets(self, dataid):
        self.datadict = {}
        data_strs = exp_number_to_str(dataid.data_nrs)
        for i in range(len(data_strs)):
            path = self.data_folder + dataid.filenamestart + data_strs[i] + dataid.ext
            time, counts = np.loadtxt(path, unpack=True)
            if self.measured_chopper_freq:
                path = self.data_folder + dataid.chopperfilename + data_strs[i] + '.txt'
                chopper_freq = np.loadtxt(path, usecols=1, unpack=True)
                chopper_freq[chopper_freq > 3*self.chopper_freq] /= 2
                chopper_freq /= 2
                chopper_freq = np.average(chopper_freq)
            else:
                chopper_freq = self.chopper_freq
            self.datadict[dataid.data_nrs[i]] = Dataset(time, counts, dataid.data_positions[i], dataid.data_vext[i],
                                                              chopper_freq)  

    def plot_raw_data(self):
        for nr in self.datadict.keys():
            plt.plot(self.datadict[nr].time, self.datadict[nr].counts, label = str(nr))
        plt.axis([0.55, 0.7, 10000, 25000])
        plt.legend()
        plt.xlabel('Time (ms)')
        plt.ylabel('Counts')
        plt.show()
        plt.close()

            
    def calc_v(self): #calculate speed of molecules using different UTI positions
        """
        
        """
        self.peaktimelist = []
        poslist = []
        
        for peaknr in range(int(np.ceil(self.len_dataset*4)/2)):
            peaktimes = []
            first = True
            for dataset in self.datadict.values():
                chopper_correction = dataset.chopper_freq / list(self.datadict.values())[0].chopper_freq 
#                chopper_correction = 1 #because this is not the right way to correct for chopper, fix later
                if peaknr==0:#not the nicest way of programming this probably
                    poslist.append(dataset.position)     
                if first:
                    peaktime0 = dataset.peaklist[peaknr].tmax
                    peaktimes.append(0)
                    first=False
                else:
                    peaktimes.append((dataset.peaklist[peaknr].tmax-peaktime0) * chopper_correction)
            print (peaktimes)
            self.peaktimelist.append(np.array(peaktimes))
            
#        print (poslist, self.peaktimelist)

        for peaktimes in self.peaktimelist:
            v, p1 = np.polyfit(peaktimes, np.array(poslist), deg=1) #datadict.keys has to be changed
            v = -v #convention
            plt.plot(peaktimes, np.array(poslist),marker='o', label=str(int(v))+' m/s') #datadict.keys
        plt.legend()
        plt.show()
        plt.close()
        
    def calc_Tion(self):#calculates Tion using different extraction potentials
        """
        Assumptions: 1+ ion charge, perfect expansion, ionisation at center of filament
        Also, fitting a straight line through something that was converted from 
        something not straight, favors certain datapoints over others (low/high values)
        """
        print ('not finished')
    
        self.peaktimelist = [] #list of Ttot
        vlist = [] #list of Vext 
        
        for peaknr in range(int(np.ceil(self.len_dataset*4)/2)):
            peaktimes = []
            for dataset in self.datadict.values():
                chopper_correction = dataset.chopper_freq / list(self.datadict.values())[0].chopper_freq 
#                peaktimes.append(dataset.peaklist[peaknr].time[dataset.peaklist[peaknr].approxpos] * chopper_correction)
                peaktimes.append(dataset.peaklist[peaknr].tmax * chopper_correction)
                vlist.append(dataset.vext)  
            
            self.peaktimelist.append(np.array(peaktimes))

        Emol = self.gas_mass / self.avg_mass * 2.5 * 8.314 * self.nozzle_temp #approximate kinetic energy of the molecules before ionisation (in J/mol)    
        Eion = np.array(vlist)*96*1000 #eV to J/mol
        
        inv_ion_speed = np.sqrt(self.gas_mass / (2*(Emol+Eion))) #kg/mol / (kgm^2s-2/mol) = s^2/m^2 -> inverse speed in s/m or ms/mm
  
        acc_times = ((np.sqrt(2*Emol/self.gas_mass) + np.sqrt((2*Emol/self.gas_mass)
                      +2*np.array(vlist)*96483/self.gas_mass))/(np.array(vlist)*96483/0.015/self.gas_mass)) #acceleration time in seconds
        
        acc_times *= 1000
        
        expected_ion_times = 0.21*inv_ion_speed *1000
        
#        print (acc_times)
        inv_ion_speed = vlist #for testing only
        
        predicted_times = predict_tion(vlist, nozzle_temp=self.nozzle_temp, gas_mass=self.gas_mass, avg_mass=self.avg_mass)

        for peaktimes in self.peaktimelist:
#            L, p1 = np.polyfit(inv_ion_speed,peaktimes-acc_times, deg=1) #datadict.keys has to be changed
#            fit = L*inv_ion_speed + p1
            plt.plot(inv_ion_speed,peaktimes,marker='o', label='original')
#            plt.plot(inv_ion_speed,peaktimes-acc_times,marker='o', label='corrected for acc') #datadict.keys
#            plt.plot(inv_ion_speed,fit, label=str(int(L))+' mm ion path length')
            plt.plot(inv_ion_speed,acc_times+expected_ion_times+peaktimes[np.argmin(inv_ion_speed)]-acc_times[np.argmin(inv_ion_speed)]-expected_ion_times[np.argmin(inv_ion_speed)],label='Expected (for original)')
            plt.plot(inv_ion_speed,predicted_times+peaktimes[np.argmin(inv_ion_speed)]-predicted_times[np.argmin(inv_ion_speed)],label='Predicted with acc and dec')
         
        plt.xlabel('inverse ion velocity (ms/mm)')
        plt.ylabel('total flight time (ms)')
        plt.legend()
        plt.show()
        plt.close()     
        
    def calc_Tslit(self): #calculates whether the distance between the two slits is the same on both sides
        """
        This only works when it triggers on the same slit every time 
        (so number of cycles = x.875 and not x.375)
        """
        t1list = []
        t2list = []
        for dataset in self.datadict.values():
            peaktimes = []
            for peaknr in range(int(np.ceil(self.len_dataset*4)/2)):
                peaktimes.append(dataset.peaklist[peaknr].tmax)
            peaktimes = np.array(peaktimes)
#            print ('peaktimes ', peaktimes)
            t1 = peaktimes[1::2]-peaktimes[0:-1:2]  
            t2 = peaktimes[2:-1:2]-peaktimes[1:-2:2]    
#            print ('t1 = ', t1)
#            print ('t2 = ', t2)  
            t1list.append(t1)
            t2list.append(t2)
        t1list = np.array(t1list)
        t2list = np.array(t2list)
        print ('t1 avg ', np.average(t1list), ' std ',np.std(t1list))
        print ('t2 avg ', np.average(t2list), ' std ',np.std(t2list))
        dt = np.absolute(np.average(t1list)-np.average(t2list))
        print ('dt = ', dt)
            
class Dataset:
    """Eigenschappen om later toe te voegen: piekposities? of losse
    class voor pieken, met alle fitparameters erin"""
    def __init__(self, time, counts, position, vext, chopper_freq, firstpeak='small'):
        self.time = time
        self.counts = counts
        self.position = position
        self.vext = vext
        self.chopper_freq = chopper_freq
        self.firstpeak = firstpeak
        
    def makepeaks(self, experiment):
        npeaks = experiment.len_dataset * 4
        window = int(len(self.time) / npeaks) #window size in data for each peak
        if self.firstpeak == 'small':
            peaknrs = np.arange(0,npeaks, 2).astype(int) #to select the small peaks only
        else:
            peaknrs = np.arange(1, npeaks+1, 2).astype(int)
            
        self.peaklist = [] 
        for peaknr in peaknrs:
            starttime = peaknr * window #index where the window starts  
            peaktime = self.time[starttime+1:starttime+window] #+1 because division by 0 later
            peakcounts = self.counts[starttime+1:starttime+window]
            self.peaklist.append(Peak(self.time[starttime], peaktime, peakcounts, peaknr))         
     
    def fitpeaks(self, experiment):
        for peak in self.peaklist:
            peak.fitpeakshape(experiment, self)
        
class Peak:
    def __init__(self, starttime, time, counts, peaknr):
        self.starttime = starttime 
        self.approxpos = np.argmax(counts)
        self.time = time
        self.counts = counts
        self.peaknr = peaknr
        
    def fitpeakshape(self, experiment, dataset):
        """
        calculates tmax = most probable time in distribution corrected for density 
        sensitive detector
        """
        A = 1e-12 * np.min(self.counts) #this factor is just a guess
        B = np.min(self.counts) #base line
        L = 567 - dataset.position #path length
        ts = self.starttime + 0.00999999 #start time = time at the start of the window, accounts for any shifts in the data
        alpha = np.sqrt(2 * 8.314 * experiment.nozzle_temp / experiment.gas_mass) #width of the peak (but also position)
        t0 = experiment.t0_guess #fit is very sensitive to t0. t0_guess has to be found by trial and error (but is usually close to peak)
        t02 = (2*(L/alpha)**2*(np.argmax(self.counts)-ts))/(2*(L/alpha)**2-4*(np.argmax(self.counts)-ts)**2)
        print ('t02', t02)

#        print ('tmax ', tmax)
#        print ('ts ', ts)
#        print ('t0 ', t0)        
        parameters = np.array([A, B, L, t0, ts, alpha, dataset.chopper_freq])

        factors = np.array([1e10, 10, 1.01, 2, 2, 2, 1.001])
#        factors = np.full(7, 1e50)
        parmin = parameters/factors
        parmin[4] -= 1
        parmax = parameters*factors
        parmax[4] += 1
        
        self.fitparams, self.fitsd = curve_fit(fitfunction_convolution, 
                                               self.time, self.counts, 
                                               p0=parameters, bounds=(parmin, parmax), maxfev=3000) #, bounds=(parmin, parmax)

        self.t0 = self.fitparams[3]
        self.ts = self.fitparams[4]
        tfit = np.arange(np.min(self.time[1:]), np.max(self.time), 0.00001) #for more precise peak position than with self.time
        yfit = fitfunction_convolution(tfit, *self.fitparams, corr_dens=True)
        self.tmax = tfit[np.argmax(yfit)] 
        yfit -= np.min(yfit)
        self.tavg = np.sum(tfit*yfit)/np.sum(yfit)
        print ('tavg ',self.tavg)
        print ('ts', self.ts)
        print ('t0', self.t0)
        
        if plotpeakshape:
            yfit = fitfunction_convolution(self.time, *self.fitparams)
            yguess = fitfunction_convolution(self.time, *parameters)
            yfit_corr = fitfunction_convolution(self.time[1:], *self.fitparams, corr_dens=True)
            print (np.max(yfit_corr))
#            yfit_corr = yfit_corr / np.max(yfit_corr) * np.max(yfit)
            print (self.fitparams)
            plt.plot(self.time, yguess, label='Guess')
            plt.plot(self.time, self.counts)              
            plt.plot(self.time, yfit, label='Fit')
            plt.plot(self.time[1:], yfit_corr, label='Fit, corr for dens')
            plt.scatter(self.tmax,np.max(self.counts))
            plt.axis([None,None,np.min(self.counts),np.max(self.counts)])
            plt.xlabel('Time (ms)')
            plt.legend()
            plt.show()
            plt.close()
        else:
            yfit = fitfunction_convolution(self.time, *self.fitparams)
            yfit -= np.min(yfit)
            yfit /= np.max(yfit)
            plt.plot(self.time, yfit, label=str(dataset.vext))
#        
        
############### Functions ####################################################
        
def exp_number_to_str(numbers, filenamestart=''):
    """numbers: array or list with numbers"""
    strs = []
    for i in range(len(numbers)):
        name = filenamestart + str(int(numbers[i]/10)) + str(int(numbers[i]%10))
        strs.append(name)
    return strs
        

def chopperfunction(x_values, r_beam, w_slit):
    """x_values: if you want n amplitudes, length x_values is ceil(1/2n+1) and numbers range from -1/2w_slit to R"""
    
    edge = x_values >= r_beam - w_slit
    center = x_values < r_beam - w_slit
    
    y_values = np.zeros(len(x_values))
    
    def integral(r_beam, x):
        return x*np.sqrt(r_beam**2-x**2) + r_beam**2*np.arcsin(x/r_beam)
    
    y_values[edge] = integral(r_beam, r_beam) - integral(r_beam, x_values[edge])
    y_values[center] = (integral(r_beam, x_values[center]+w_slit) 
                        - integral(r_beam, x_values[center]))
    
    return y_values
       
def get_chopper_amplitudes(n, chopper_freq, r_beam=0.23, w_slit=0.85): #r and w in mm
    """ Returns only odd length array. if n is even, returns array of length n-1"""
#    xt_conversion = 24 * 250 / (w_slit * chopper_freq) #in microseconds, old conversion that is not correct
    xt_conversion = 1/(2*np.pi*55.992*chopper_freq) * 1e3 #1e3 to convert to milliseconds

    
    x_values = np.linspace(-w_slit/2, r_beam, num=int(np.ceil(n/2+1)))
    
    amplitudes = chopperfunction(x_values, r_beam, w_slit)
    amplitudes = np.concatenate((np.flip(amplitudes)[:-1],amplitudes))
    amplitudes = amplitudes[1:-1]
    amplitudes /= np.pi * r_beam**2 #to normalize
    amplitudes[np.isnan(amplitudes)] = 1 #very ugly, but it works..

    
    t_values = (x_values+w_slit/2) * xt_conversion
    t_values = np.concatenate(((-1*np.flip(t_values)[:-1]),t_values))
    t_values = t_values[1:-1]
    
    return t_values, amplitudes

def fitfunction(t, A, B, L, t0, ts, alpha, chopper_freq):
    return ( A * ((L/(t-ts))**4 
              * np.exp(-((L/(t-ts)-L/(t0))/alpha)**2)) + B)

def fitfunction_convolution(t, A, B, L, t0, ts, alpha, chopper_freq,corr_dens=False):
    """
    A = amplitude
    B = background
    L = flight path length
    t0 = peak maximum position
    ts = time shift of the function
    alpha = width of the peak
    """
    n = 15
#    print (chopper_freq)
    t_shift, amp = get_chopper_amplitudes(n, chopper_freq)
#    print (t_shift)
    
    y = np.zeros(len(t))
    for i in range(n):
        if corr_dens:
            y += (amp[i] * A * ((L/(t-ts-t_shift[i]))**4 * np.exp(-((L/(t-ts-t_shift[i])-L/t0)/alpha)**2)) / t) #
        else:
            y += (amp[i] * A * ((L/(t-ts-t_shift[i]))**4 * np.exp(-((L/(t-ts-t_shift[i])-L/t0)/alpha)**2))) #
    y /= np.sum(amp)  
    y += B
    
    return y 


def testplot_fitfunctions():
    t = np.arange(0.001, 10, 0.00001)
    A = 1e-13
    B = 0
    L = 563
    v0 = 2000

    ts = 0.0999999
    alpha = 1000 * 1
    
    tmaxlist = []
    tmaxlist_corr = []
    L_list = [0.000001, 500, 1000, 5000, 10000]
    
    for L in L_list:
        t0 = L/v0
        y = fitfunction(t, A, B, L, t0, ts, alpha, 1)
        y2 = y/t
        y /= np.max(y)
        y2 /= np.max(y2)
        tmax = t[np.argmax(y)]
        tmax_corr = t[np.argmax(y2)]
        tmaxlist.append(tmax)
        tmaxlist_corr.append(tmax_corr)
        plt.plot(t, y, label = 'density detector, tmax = '+str(np.round(tmax,4)))
        plt.plot(t, y2, label = 'corr, tmax = '+str(np.round(tmax_corr,4)))
    plt.title('t0 = '+str(t0))
    plt.axis([None,None,0,1])
    plt.legend()
    plt.show()
    plt.close()
    
    L=0.00000000001
    t0 = L/v0
    
    for ts in [0.0099999, 0.049999999, 0.099999999, 0.4999999999]:
        print ('ts', ts)
        tmax = t[np.argmax(fitfunction(t, A, B, L, t0, ts, alpha, 1))]
        print ('tmax at L->0', tmax)
        
    
    #Check linearity of peak position vs L 
    slope, y0 = np.polyfit(np.array(L_list[1:]),np.array(tmaxlist[1:]),1)
    fit = y0 + slope*np.array(L_list)
    slope, y0 = np.polyfit(np.array(L_list[1:]),np.array(tmaxlist_corr[1:]),1)
    fit_corr = y0 + slope*np.array(L_list)
    
    plt.plot(L_list, tmaxlist,label='density detector')
    plt.plot(L_list, fit, label = 'fit (dens) ')
    plt.plot(L_list, tmaxlist_corr, label='corrected')
    plt.plot(L_list, fit_corr, label = 'fit (corr)')
    plt.legend()
    plt.show()
    plt.close()
    
    

    for L in [500, 1000, 2000, 10000]:
        t0 = L/v0    
        v = L/(t-ts)
        yv = fitfunction(t, A, B, L, t0, ts, alpha, 1) * 1/v**2
        yv2 = yv*v
        yv /= np.max(yv)
        yv2 /= np.max(yv2)
        plt.plot(v, yv,label='density, vmax = '+str(v[np.argmax(yv)]))
        plt.plot(v, yv2,label='corr, vmax = '+str(v[np.argmax(yv2)])+', vavg = '+str(np.sum(v*yv2)/np.sum(yv2)))
    plt.title('v0 = ' + str(v0))
    plt.legend()
    plt.axis([0,20000,0,1])
    plt.show()
    plt.close()
    

    
def testplot_chopper(exp):
    t, y = get_chopper_amplitudes(25, exp.chopper_freq)
    plt.plot(t*1000, y)
    plt.xlabel('Time (us)')
    plt.show()
    plt.close()

def testplot_tion(vlist, vfoc=20, nozzle_temp=300, gas_mass=0.004, avg_mass=0.004, Lion=0.21):

    Emol = 2.5 * 8.314 * nozzle_temp #approximate kinetic energy of the molecules before ionisation (in J/mol)    
    Eion = np.array(vlist)*96*1000 #eV to J/mol, kinetic energy ions gain during acceleration
    
    inv_ion_speed = np.sqrt(gas_mass / (2*(Emol+Eion))) #kg/mol / (kgm^2s-2/mol) = s^2/m^2 -> inverse speed in s/m or ms/mm
  
    acc_times = ((np.sqrt(Emol/2/gas_mass) + np.sqrt((Emol/2/gas_mass)
                  +2*np.array(vlist)*2.4088e7))/(np.array(vlist)/0.015*2.4088e7)) #acceleration time in seconds    
    acc_times *= 1000 #in ms
    
    ion_times = Lion*inv_ion_speed
    ion_times *= 1000 #in ms
    
    plt.plot(vlist, acc_times+ion_times)
    plt.show()
    plt.close()
    
def predict_tion(Vlist, Vfoc=20, nozzle_temp=300, gas_mass=0.004, avg_mass=0.004, Lion=0.21, Lfil=0.005, ms=True):
    Emol = 2.5*8.314*nozzle_temp
    vmol = np.sqrt(2*Emol/avg_mass)  
    
    acc_times = ((-vmol + np.sqrt(vmol**2+2*(np.array(Vlist)+Vfoc)*96483/gas_mass))
                /((np.array(Vlist)+Vfoc)*96483/Lfil/gas_mass))
    
    v0 = np.sqrt(2*(gas_mass/avg_mass*Emol+(np.array(Vlist)+Vfoc)*96*1000)/gas_mass)
    
    dec_times = ((-v0 + np.sqrt(v0**2+2*(1000)*96483/gas_mass))
                /((1000)*96483/Lion/gas_mass))
    
    if ms:
        acc_times *= 1000
        dec_times *= 1000
    
    print ('acc ',acc_times)
    print ('dec ',dec_times)
    
    return acc_times+dec_times

main()