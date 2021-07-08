# -*- coding: utf-8 -*-
"""
This script uses the first peak in each TOF dataset (assumes it is a small peak)
to calculate the flight time of the molecules

https://mail.python.org/pipermail/scipy-user/2013-April/034403.html

@author: Charlotte
"""
import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit, least_squares

plotpeakshape = True

def main():
########### Settings ########################################################     
    
    data_folder = 'P:/DATA/2019/04 Apr/190430/TOF/'
#    data_folder = 'P:/DATA/2019/07 Jul/190725/TOF/'
#    data_folder = 'P:/DATA/2019/08 Aug/190808/TOF/'
#    data_folder = 'P:/DATA/2019/08 Aug/190812/TOF/'
#    data_folder = 'P:/DATA/2019/05 May/190502/TOF/'
    data_folder = 'P:/DATA/2020/01 Jan/200113/TOF/'
    data_folder = 'P:/DATA/2020/01 Jan/200128/TOF/'
#    data_folder = 'C:/Users/Charlotte/surfdrive/DATA/TOF testdata/'   
#    data_folder = 'P:/DATA/TOF testdata/New WinRAR archive/170612/tof/'

    nr_of_datasets = 7
    start_dataset = 0
    data_nrs = np.linspace(start_dataset, nr_of_datasets+start_dataset-1, num=nr_of_datasets) #numbers of datasets

    data_positions = np.full(nr_of_datasets,46)   
    data_positions = [0, 46, 35, 25, 15, 5, 0] #positions on moving stage (mm), corresponding to number of dataset above
    
    measured_chopper_freq = False
    chopper_freq = 250 #Hz
    gas_mass = 0.004 #kg/mol
    avg_mass = (4*0.004+1*0.04+2.9*0.002)/7.9 #calculate using beam composition
    avg_mass = (1*0.004+2.5*40)/3.5
    nozzle_temp = 300 #K


    data_vext = np.full(nr_of_datasets, 15)



############# Script ########################################################  
    dataid = Dataidentifier(data_nrs, data_positions, data_vext)
     
    exp = Experiment(data_folder, dataid, chopper_freq, gas_mass, avg_mass, nozzle_temp, measured_chopper_freq)
    exp.plot_raw_data()
  
    exp.fitparams = exp.fitpeaks()
   
    for nr, dataset in zip(range(len(exp.datadict)),exp.datadict.values()):
        print (nr)
        dataset.plotfit(exp, nr)
    
    exp.plot_v_E_distr()
    


    
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
        
        #Guessed parameters of the time/velocity distribution
        self.L_guess = 567
        
        self.v_guess = np.sqrt(2*2.5*8.314*self.nozzle_temp/self.avg_mass) #for perfect expansion
        self.t_window = self.L_guess/self.v_guess + 1000/8/self.chopper_freq #approximate peak position + 1/8 of chopper cycle. in ms
        self.index_window = int(self.t_window * 2000) #assuming 500 ns bins
        
        self.plot_raw_data()
        
        self.ts_guess, L_tmax, tmax_guess = self.guess_ts_tmax()
        self.alpha_guess = np.sqrt(2 * 8.314 * self.nozzle_temp /self.avg_mass)*0.3 #width of the peak (but also position). physical meaning: velocity?
        self.t0_guess = (2*(L_tmax/self.alpha_guess)**2*(tmax_guess-self.ts_guess)
                       /(2*(L_tmax/self.alpha_guess)**2-4*(tmax_guess-self.ts_guess)**2))
        self.v0_guess = L_tmax / self.t0_guess
        self.A_guess = 1E-10
        self.B_guess = np.min(list(self.datadict.values())[0].counts)*1.5
        
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
                         dataid.data_positions[i], dataid.data_vext[i], chopper_freq)  

    def plot_raw_data(self):
        for nr in self.datadict.keys():
            plt.plot(self.datadict[nr].time[:self.index_window], 
                     self.datadict[nr].counts[:self.index_window], label = str(nr))
        plt.axis([None, None, None, None])
        plt.legend()
        plt.xlabel('Time (ms)')
        plt.ylabel('Counts')
        plt.show()
        plt.close()
        
    def guess_ts_tmax(self, plot=False):
        """
        Guesses ts based on a linear fit through tmax (gauss fit) vs L
        Returns ts_guess and L and tmax for calculating t0 and v0
        """   
        positions = []
        tmax = []
        
        def fit_gauss(dataset):
            x = dataset.time[:self.index_window]
            y = dataset.counts[:self.index_window]
            thresh = (np.min(y)+np.max(y))/2# select only top half of the peak
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
            
            return x[np.argmax(fit)]
        
        for dataset in self.datadict.values():
            positions.append(self.L_guess - dataset.position)
            tmax.append(fit_gauss(dataset))
        
        positions = np.array(positions)
        tmax = np.array(tmax)
        

        slope, ts_guess = np.polyfit(positions, tmax, 1)
        
        plt.plot(positions, tmax)
        plt.plot(positions, ts_guess+positions*slope, label='fit, ts = '+str(ts_guess))
        plt.ylabel('time (ms)')
        plt.xlabel('Position (mm)')
        plt.legend()
        plt.show()
        plt.close()
        
        return ts_guess, positions[0], tmax[0]

    def fitpeaks(self):
   
        def err_multiple_leastsq(parameters):
#            print (parameters)
            residual_list=np.array([])
            for nr, dataset in zip(range(len(self.datadict)),self.datadict.values()):
                parameters2 = parameters[[int(nr), int(nr)+len(self.datadict), -4,-3,-2,-1]]
                residual_list = np.concatenate((residual_list, err_leastsq(parameters2, dataset.time[:self.index_window], dataset.counts[:self.index_window], dataset.chopper_freq, dataset.position)))
                
#            print(residual_list)
            return residual_list[1:]
        
        factors = [50 for i in range(len(self.datadict))] + [5 for i in range(len(self.datadict))]
        factors += [1.01, 2, 10, 5]
        factors = np.array(factors)        
        
        bounds = (self.params/factors, self.params*factors)
        result = least_squares(err_multiple_leastsq, self.params, bounds=bounds, verbose=1) #, x_scale=1/np.array([1e-10, 1e2, 500, 0.3, 0.01, 1e3])
        
        print ('success',result.success)
        print (result.message)
        print ('parameters', self.params)
        print ('x', result.x)
        print ('ts', result.x[-2])
        return result.x
        
    def plot_v_E_distr(self):
        """
        This should be done without a baseline and time shift.
        Not sure of this is mathematically correct at this moment.
        Does correct for density sensitive detector
        """
        vmax = 4000
        v = np.arange(1,vmax,1)
        t = self.fitparams[-4]/v - self.fitparams[-2]
        y = fitfunction_convolution(t, 100, 0, self.fitparams[-4], self.fitparams[-3],0,self.fitparams[-1] ,self.chopper_freq,corr_dens=True,use_v0=True) 
        yv = y / v**2
        plt.plot(v,yv/np.sum(yv))
        plt.title('Velocity distribution, v_avg = '+str(np.round(np.sum(v*yv)/np.sum(yv)))+' m/s')
        plt.ylabel('Probability density')
        plt.xlabel('v (m/s)')
        plt.axis([0,vmax,0,None])
        plt.show()
        plt.close()
        
        E = 0.5 * self.gas_mass * v**2 / 1000 #kJ/mol
        yE = 1 / (self.gas_mass * v) * yv
        
        plt.plot(E, yE/np.sum(yE))
        plt.title('Energy distribution, E_avg = '+str(np.round(np.sum(E*yE)/np.sum(yE),3))+' kJ/mol')
        plt.ylabel('Probability density')
        plt.xlabel('E (kJ/mol)')
        plt.axis([0,np.max(E),0,None])
        plt.show()
        plt.close()
        
            
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

     
    def plotfit(self, experiment, nr):
        guess = fitfunction_convolution(self.time[:experiment.index_window],
                *experiment.params[[int(nr), int(nr)+len(experiment.datadict), -4,-3,-2,-1]], self.chopper_freq, pos = self.position,use_v0=True)
        
        fit = fitfunction_convolution(self.time[:experiment.index_window],
              *experiment.fitparams[[int(nr), int(nr)+len(experiment.datadict), -4,-3,-2,-1]], self.chopper_freq, pos = self.position,use_v0=True)
        
        plt.plot(self.time[:experiment.index_window],self.counts[:experiment.index_window], label='Data')
        plt.plot(self.time[:experiment.index_window],guess, c='g', label='Guess')
        plt.plot(self.time[:experiment.index_window],fit, c='r',label='Fit')
        plt.axis([0,None,0,np.max(self.counts[:experiment.index_window])])
        plt.title('pos = '+str(self.position)+' mm, max_fit = '+str(self.time[np.argmax(fit)]))
        plt.xlabel('Time (ms)')
        plt.ylabel('Counts')
        plt.legend()
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

"""OLD THINGS THAT ARE NOT USED ANYMORE""""""

def err(parameters, t, counts, chopper_freq, pos):
#    not used anymore

#    print ('parameters err',parameters)
    y = fitfunction_convolution(t, *parameters, chopper_freq, pos=pos)  
#    print (y)
    return np.sum((y-counts)**2)
    
    
FROM WITHIN EXPERIMENT.FITPEAKS:
    
        def err_multiple(parameters): #not used anymore
            som = 0
            for dataset in self.datadict.values():
#                print ('parameters err_multiple \n', parameters)
                som += err(parameters, dataset.time[:self.index_window],
                           dataset.counts[:self.index_window], dataset.chopper_freq, dataset.position)
#            print (som)
            return som
        
#        bounds = []
#        for i in range(len(self.params)):
#            bounds.append((self.params[i]/factors[i], self.params[i]*factors[i]))
#        result = minimize(err_multiple, self.params, args=(), bounds=bounds) #, method='L-BFGS-B', options={'ftol':1E-20, 'maxiter':10000000, 'maxfun':10000000, 'factr':10}        

END OF PART FROM EXPERIMENT.FITPEAKS


    def calc_v(self): #calculate speed of molecules using different UTI positions

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
        
"""