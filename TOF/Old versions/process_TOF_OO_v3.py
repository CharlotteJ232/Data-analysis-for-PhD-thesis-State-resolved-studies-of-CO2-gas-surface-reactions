# -*- coding: utf-8 -*-
"""
This script uses the first peak in each TOF dataset (assumes it is a small peak)
to calculate the flight time of the molecules

https://mail.python.org/pipermail/scipy-user/2013-April/034403.html

New in v3: now also works if the beam is so slow that a part of the big peak is 
visible at the start of the dataset

@author: Charlotte
"""
import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit, least_squares
import os


########### Settings ########################################################

save_all_plots = True

directory = 'P:/Surfdrive/'
#directory = 'D:/Surfdrive/'
  
#data_folder = 'P:/DATA/2019/12 Dec/191223/TOF/'
#avg_mass = (4*0.004+1*0.04+2.9*0.002)/7.9

data_folder = directory+'DATA/2020/01 Jan/200128/TOF/'
avg_mass = (0.1*0.004+2.5*0.040)/3.5

data_folder = directory+'DATA/2020/02 Feb/200204/TOF/'
avg_mass = 0.004

#data_folder = directory+'DATA/2020/02 Feb/200207/TOF/'
#avg_mass = (0.5*0.004+3*0.040)/3.5
#
#data_folder = directory+'DATA/2020/02 Feb/200211/TOF/'
#avg_mass = (1*0.004+8*0.020)/9
#
#data_folder = directory+'DATA/2020/02 Feb/200218/TOF/'
#avg_mass = (9*0.004+1*0.040)/10
##
#data_folder = directory+'DATA/2020/02 Feb/200220/TOF/'
#avg_mass = (3*0.004+2*0.040)/10
#
#data_folder = directory+'DATA/2020/03 Mar/200306/TOF/'
#avg_mass = (19*0.004+1*0.040)/20
##
#data_folder = directory+'DATA/2020/03 Mar/200310/TOF/'
#avg_mass = (0.5*0.004+3*0.040)/3.5
##
#data_folder = directory+'DATA/2020/03 Mar/200311/TOF/'
#avg_mass = (19*0.004+1*0.040)/20
#
#data_folder = directory+'DATA/2020/05 May/200511/TOF/'
#avg_mass = (8*0.004+2*0.040)/10
#
#data_folder = directory+'DATA/2020/05 May/200519/TOF/'
#avg_mass = (19*0.004+0.5*0.040)/19.5

data_folder = directory+'DATA/2020/07 Jul/200724/TOF/'
avg_mass = 0.004

data_folder = directory+'DATA/2020/07 Jul/200727/TOF/'
avg_mass = 0.004

data_folder = directory+'DATA/2020/08 Aug/200813/TOF/'
avg_mass = 0.004

#data_folder = directory+'DATA/2020/08 Aug/200813_2/TOF/'
#avg_mass = 0.004

data_folder = directory+'DATA/2020/08 Aug/200825/TOF/'
avg_mass = 0.004

data_folder = directory+'DATA/2020/09 Sep/200904/TOF/'
avg_mass = 0.004

start_dataset = 0
nr_of_datasets = 7

gas_mass = 0.004 #kg/mol
#if gas_mass == 0.004:
#    start_dataset=0
#if gas_mass == 0.04:
#    start_dataset=7



data_nrs = np.linspace(start_dataset, nr_of_datasets+start_dataset-1, num=nr_of_datasets) #numbers of datasets

data_positions = np.full(nr_of_datasets,46)   
data_positions = [0, 46, 35, 25, 15, 5, 0] #positions on moving stage (mm), corresponding to number of dataset above

measured_chopper_freq = False
chopper_freq = 250 #Hz
nozzle_temp = 600 #K

savefolder = data_folder+'Images/'+str(gas_mass * 1000)+'/'

data_vext = np.full(nr_of_datasets, 15) #at this moment, this does nothing
increase_fit_bounds = 1 #increases the fitparameter bound range by this factor. If you use this, check if the path length is not fitted to an incorrect value

def main():

############# Script ########################################################  
    dataid = Dataidentifier(data_nrs, data_positions, data_vext)
     
    exp = Experiment(data_folder, dataid, chopper_freq, gas_mass, avg_mass, nozzle_temp, measured_chopper_freq)
    exp.plot_raw_data(save=save_all_plots)
  
    exp.fitparams = exp.fitpeaks()
   
    for nr, dataset in zip(range(len(exp.datadict)),exp.datadict.values()):
        print (nr)
        dataset.plotfit(exp, nr, save=save_all_plots,isolate_peak=True)
    
    exp.plot_v_E_distr(save=save_all_plots)
    


    
################# Classes ####################################################


class Dataidentifier:
    def __init__(self, data_nrs, data_positions, data_vext, filenamestart='TOF', 
                 chopperfilename = 'Chopper', ext='.asc'):
        self.data_nrs = data_nrs
        self.data_positions = data_positions
        self.data_vext = data_vext
        self.filenamestart = filenamestart
        self.chopperfilename = chopperfilename
        self.ext = ext


class Experiment:
    def __init__(self, data_folder, dataid, chopper_freq, 
                 gas_mass, avg_mass, nozzle_temp, measured_chopper_freq):
        self.chopper_freq = chopper_freq
        self.gas_mass = gas_mass
        self.avg_mass = avg_mass
        self.nozzle_temp = nozzle_temp
        self.data_folder = data_folder
        self.measured_chopper_freq = measured_chopper_freq
        
        self.read_datasets(dataid)
        
#        self.plot_raw_data()
        
        #Guessed parameters of the time/velocity distribution
        self.L_guess = 567
        
        self.E_guess = 2.5*8.314*self.nozzle_temp #for perfect expansion
        self.v_guess = np.sqrt(2*self.E_guess/self.avg_mass) #for perfect expansion
        self.t_guess = self.L_guess/self.v_guess   #approximate peak position  
        self.t_window = 1000/12/self.chopper_freq     #1/12 of chopper cycle. in ms. this will be added to the left and right of the peak to create the data window
        self.index_guess = int(self.t_guess * 2000) #assuming 500 ns bins
        self.index_window = int(self.t_window * 2000) #assuming 500 ns bins
        if self.index_window > self.index_guess:
            self.index_window = self.index_guess
        
        self.ts_guess, L_tmax, tmax_guess = self.guess_ts_tmax(plot=True) #ts = t_shift, which is all delay times in the data combined (electronic, mass spectrometer, etc)
        self.alpha_guess = np.sqrt(2 * 8.314 * self.nozzle_temp /self.avg_mass)*0.1 #width of the peak (but also position). physical meaning: velocity?
        self.t0_guess = (2*(L_tmax/self.alpha_guess)**2*(tmax_guess-self.ts_guess) #not sure of the physical meaning, derived this from the equation.
                       /(2*(L_tmax/self.alpha_guess)**2-4*(tmax_guess-self.ts_guess)**2))
        self.v0_guess = L_tmax / self.t0_guess #stream velocity, whatever that means. this is used instead of t0, because it is independent of path length and can thus be fitted simultaneously for all datasets
        self.A_guess = 1E-10 #amplitude, random factor
        self.B_guess = np.min(list(self.datadict.values())[0].counts)*1.5 #background level
        
        self.params = [self.A_guess for i in range(len(self.datadict))] 
        self.params += [self.B_guess for i in range(len(self.datadict))]
        self.params += [self.L_guess, self.v0_guess, self.ts_guess, self.alpha_guess]
        self.params = np.array(self.params)
        

        
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
            self.datadict[dataid.data_nrs[i]] = Dataset(time, counts, 
                         dataid.data_positions[i], dataid.data_vext[i], dataid.data_nrs[i], chopper_freq)  

    def plot_raw_data(self,isolate_peak = False, save=False):
        for nr in self.datadict.keys(): #plot all datasets in one figure
            if isolate_peak: #plot only the peak we are interesed in
                plt.plot(self.datadict[nr].time[self.index_guess-self.index_window:self.index_guess+self.index_window], 
                     self.datadict[nr].counts[self.index_guess-self.index_window:self.index_guess+self.index_window], label = str(nr))
            else: #plot entire dataset
                plt.plot(self.datadict[nr].time, 
                     self.datadict[nr].counts, label = str(nr))
        plt.axis([0,np.max(self.datadict[nr].time), 0, None])
        plt.legend()
        plt.title('Raw data')
        plt.xlabel('Time (ms)')
        plt.ylabel('Counts')
        if save:
            if not os.path.exists(savefolder):
                os.makedirs(savefolder)
            plt.savefig(savefolder+'raw.png', dpi=1000)
        plt.show()
        plt.close()
        
    def guess_ts_tmax(self, plot=False, save=False):
        """
        Guesses ts based on a linear fit through tmax (determined by gauss fit) 
        vs L.
        Returns ts_guess and L and tmax for calculating t0 and v0
        
        if plot: plots gaussian fit for each dataset, and also tmax as a function
        of position (including fit to extract ts)
        """   
        positions = []
        tmax = []
        
        def fit_gauss(dataset):

            x = dataset.time[self.index_guess-self.index_window:self.index_guess+self.index_window] #select only the peak within the window
            y = dataset.counts[self.index_guess-self.index_window:self.index_guess+self.index_window]
            thresh = (1/4*np.min(y)+3/4*np.max(y))# select only top quarter of the peak

            index = y > thresh
            x = x[index]
            y = y[index]
            
            fwhm = x[-1]-x[0]
            mean = x[0] + fwhm / 2
            
            def gauss_function(x, a, x0, sigma):
                return a*np.exp(-(x-x0)**2/(2*sigma**2))
            
            popt, pcov = curve_fit(gauss_function, x, y, p0=[np.max(y), mean, fwhm/2.355])
            fit = gauss_function(x,*popt)
            
            if plot:
                plt.plot(x,y)
                plt.plot(x,fit)
                plt.title(str(popt))
                plt.show()
                plt.close()
                
                if save:
                    if not os.path.exists(savefolder):
                        os.makedirs(savefolder)
                    plt.savefig(savefolder+'peak_gauss_fit_'+str(dataset.position)+'.png', dpi=1000)
            
            return x[np.argmax(fit)]
        
        for dataset in self.datadict.values():
            positions.append(self.L_guess - dataset.position)
            tmax.append(fit_gauss(dataset))
        
        positions = np.array(positions)
        tmax = np.array(tmax)
        

        slope, ts_guess = np.polyfit(positions, tmax, 1)
        
        if ts_guess <= 0:
            #setting to 0 and to 0.1 does not work, because it will lead to division by 0 later in the program (at the convolution with the chopper hole)
            print ('Warning! Ts is negative. Will be set to 0.00001')
            ts_guess = 0.00001
        
        if plot:
            plt.scatter(positions, tmax)
            plt.plot(positions, ts_guess+positions*slope, label='fit, ts = '+str(ts_guess))
            plt.ylabel('time (ms)')
            plt.xlabel('Position (mm)')
            plt.legend()
            plt.show()
            plt.close()
            
            if save:
                if not os.path.exists(savefolder):
                    os.makedirs(savefolder)
                plt.savefig(savefolder+'ts_guess.png', dpi=1000)
        
        return ts_guess, positions[0], tmax[0]

    def fitpeaks(self):
        """
        Fits the data in the following way:
            - For each dataset, it calculates the expected values, by taking the 
            current fitparameters and filling those in in the function. The 
            function is the molecular beam functional form for a density sensitive
            detector, and convolution with the width of the chopper hole is taken
            into account by adding a number of peaks shifted slightly with respect
            to each other.
            - Then, the residuals are calculated (data - expected values) by the
            function err_leastsq()
            - These residuals are put in a list for all datasets combined by the
            function err_multiple_leastsq()
            - Residuals are then minimized by least_squares (from scipy.optimize)
            
        Notes: The fit parameters L, t0, ts, alpha are the same for all datasets. 
        fit parameters A and B vary between datasets (as these are not physical parameters)
        """
   
        def err_multiple_leastsq(parameters):
#            print (parameters)
            residual_list=np.array([])
            for nr, dataset in zip(range(len(self.datadict)),self.datadict.values()):
                parameters2 = parameters[[int(nr), int(nr)+len(self.datadict), -4,-3,-2,-1]]
                residual_list = np.concatenate((residual_list, err_leastsq(parameters2, dataset.time[self.index_guess-self.index_window:self.index_guess+self.index_window], dataset.counts[self.index_guess-self.index_window:self.index_guess+self.index_window], dataset.chopper_freq, dataset.position)))
                
#            print(residual_list)
            return residual_list[1:]
        
        
        #Make bounds for the fit parameters
        factors = [100000 for i in range(len(self.datadict))] + [5 for i in range(len(self.datadict))] #determines the bounds on the fit parameters
        factors += [1.001, 2, 5, 5]
        factors = np.array(factors) * increase_fit_bounds  #multiplied by easily accessible variable, so we can increase fit bounds if necessary
        
        upper = self.params*factors 
        upper[-2] += 0.1 #because ts is a very small number and can sometimes be negative too. A shift works better than a factor
        lower = self.params/factors
        lower[-2] -= 0.1         
      
        bounds = (lower, upper)
        
        #The fit
#        import pdb; pdb.set_trace()
        result = least_squares(err_multiple_leastsq, self.params, bounds=bounds) #verbose=1, x_scale=1/np.array([1e-10, 1e2, 500, 0.3, 0.01, 1e3])
        
#        print ('success',result.success)
#        print (result.message)
        print ('parameters', self.params) #initial guess for parameters
        print ('x', result.x) #fitted parameters
        print ('ts', result.x[-2])
        print (result.x[-4:])
        
        paramnames = np.array(['L','v0','Ts','alpha'])
        
        for_saving = np.column_stack((paramnames, self.params[-4:], result.x[-4:]))
        
        print(for_saving)
        
        np.savetxt(data_folder+'fitparams.txt',for_saving,header='Name, guess, fit',fmt='%s')
        return result.x
        
    def plot_v_E_distr(self, save=False):
        """
        This should be done without a baseline and time shift.
        Not sure of this is mathematically correct at this moment.
        Does correct for density sensitive detector
        """
        vmax = 2*self.v_guess
        v = np.arange(1,vmax,1)
        t = self.fitparams[-4]/v - self.fitparams[-2]
        y = fitfunction_convolution(t, 100, 0, self.fitparams[-4], self.fitparams[-3],0,self.fitparams[-1] ,self.chopper_freq,corr_dens=True,use_v0=True) 
        yv = y / v**2
        plt.plot(v,yv/np.sum(yv))
        plt.title('Velocity distribution, v_avg = '+str(np.round(np.sum(v*yv)/np.sum(yv)))+' m/s')
        plt.ylabel('Probability density')
        plt.xlabel('v (m/s)')
        plt.axis([0,vmax,0,None])
        
        if save:
            if not os.path.exists(savefolder):
                os.makedirs(savefolder)
            plt.savefig(savefolder+'velocity.png', dpi=1000)        
        
        plt.show()
        plt.close()
        
        E = 0.5 * self.gas_mass * v**2 / 1000 #kJ/mol
        yE = 1 / (self.gas_mass * v) * yv
        
        #Calculate FWHM
        upper = np.argwhere(yE > np.max(yE)/2)
        upper = E[upper]
        width = upper[-1]-upper[0]
        q = np.sum(E*yE)/np.sum(yE)/width
        
        
        
        plt.plot(E, yE/np.sum(yE))
        plt.title('Energy distribution, E_avg = '+str(np.round(np.sum(E*yE)/np.sum(yE),3))+' kJ/mol, Q = '+str(np.round(q[0],2)))
        plt.ylabel('Probability density')
        plt.xlabel('E (kJ/mol)')
        plt.axis([0,np.max(E),0,None])
        
        if save:
            if not os.path.exists(savefolder):
                os.makedirs(savefolder)
            plt.savefig(savefolder+'energy.png', dpi=1000)
        
        plt.show()
        plt.close()
        
        X = np.column_stack((E, yE/np.sum(yE)))
        np.savetxt(savefolder+'energy.txt',X,header='Energy (kJ/mol), probability density')
        
            
class Dataset:
    """Eigenschappen om later toe te voegen: piekposities? of losse
    class voor pieken, met alle fitparameters erin"""
    def __init__(self, time, counts, position, vext, data_nr, chopper_freq, firstpeak='small'):
        self.time = time
        self.counts = counts
        self.position = position
        self.vext = vext
        self.data_nr = data_nr
        self.chopper_freq = chopper_freq
        self.firstpeak = firstpeak    

     
    def plotfit(self, experiment, nr, save=False, isolate_peak=False):
        guess = fitfunction_convolution(self.time[experiment.index_guess-experiment.index_window:experiment.index_guess+experiment.index_window],
                *experiment.params[[int(nr), int(nr)+len(experiment.datadict), -4,-3,-2,-1]], self.chopper_freq, pos = self.position,use_v0=True)
        
        fit = fitfunction_convolution(self.time[experiment.index_guess-experiment.index_window:experiment.index_guess+experiment.index_window],
              *experiment.fitparams[[int(nr), int(nr)+len(experiment.datadict), -4,-3,-2,-1]], self.chopper_freq, pos = self.position,use_v0=True)
        
        plt.plot(self.time,
                 self.counts, label='Data')
        plt.plot(self.time[experiment.index_guess-experiment.index_window:experiment.index_guess+experiment.index_window],
                 guess, c='g', label='Guess')
        plt.plot(self.time[experiment.index_guess-experiment.index_window:experiment.index_guess+experiment.index_window],
                 fit, c='r',label='Fit')
        plt.title('pos = '+str(self.position)+' mm, max_fit = '+str(self.time[np.argmax(fit)]))
        plt.xlabel('Time (ms)')
        plt.ylabel('Counts')
        if isolate_peak:
            add = 'zoom'
            plt.axis([0,self.time[experiment.index_guess+experiment.index_window],0,np.max(self.counts[experiment.index_guess-experiment.index_window:experiment.index_guess+experiment.index_window])*1.1])
        else:
            add = 'full'
            plt.axis([0,np.max(self.time),0,np.max(self.counts[experiment.index_guess-experiment.index_window:experiment.index_guess+experiment.index_window])])
        plt.legend()
        
        if save:
            if not os.path.exists(savefolder):
                os.makedirs(savefolder)
            plt.savefig(savefolder+'fit_TOF'+str(self.data_nr)+add+'.png', dpi=1000)        
        
        plt.show()
        plt.close()
        

        
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




def err_leastsq(parameters, t, counts, chopper_freq, pos):
#    print (parameters)
    y = fitfunction_convolution(t, *parameters, chopper_freq, pos=pos, use_v0=True)  
#    print (counts,y)
    return counts-y


def fitfunction_convolution(t, A, B, L, t0, ts, alpha, chopper_freq,corr_dens=False, n=15, pos=0, use_v0=False):
    """
    A = amplitude
    B = background
    L = flight path length
    t0 = peak maximum position
    ts = time shift of the function (=delay time)
    alpha = width of the peak
    
    chopper_freq = frequency of the chopper
    corr_dens = correct for density sensitive detector (to get the actual shape 
    of the peak in the time domain, instead of what is measured)
    n = number of points for convolution
    pos = shift in L, position of mass spectrometer
    """
    t_shift, amp = get_chopper_amplitudes(n, chopper_freq) # shift/amplitude array for convolution

    L = L-pos  
    
    if use_v0:
        t0 = L/t0 #because given t0 was actually v0
  
    y = np.zeros(len(t))
    for i in range(n):
        #prevent dividing by 0
#        import pdb; pdb.set_trace()
        if np.any(t-ts-t_shift[i]==0):
            print ('dividing by 0')
            return 0

        if corr_dens:
            y += (amp[i] * A * ((L/(t-ts-t_shift[i]))**4 * np.exp(-((L/(t-ts-t_shift[i])-L/t0)/alpha)**2)) / t) #
        else:
            y += (amp[i] * A * ((L/(t-ts-t_shift[i]))**4 * np.exp(-((L/(t-ts-t_shift[i])-L/t0)/alpha)**2))) #
    y /= np.sum(amp)  #so 
    y += B
    
    return y 



main()

