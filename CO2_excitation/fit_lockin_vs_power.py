import numpy as np
from matplotlib import pyplot as plt
from matplotlib import cm 
import datetime as dt
from scipy.optimize import curve_fit, least_squares

datainfo = { 'R0':{'year':'2022', 'month':'02', 'day':'24', 'min_power':5}, #was 0.0036
            'R2':{'year':'2022', 'month':'03', 'day':'04', 'min_power':4.6},
            'R4':{'year':'2022', 'month':'03', 'day':'04', 'min_power':4.1},
            'R6':{'year':'2022', 'month':'02', 'day':'25', 'min_power':5},
            'R8':{'year':'2022', 'month':'02', 'day':'25', 'min_power':5},
            'R10':{'year':'2022', 'month':'03', 'day':'03', 'min_power':5}
            }

A_guess = { 'R0':1.0,
            'R2':0.24072303,
            'R4':0.21098072,
            'R6':0.2211775,
            'R8':0.30620438,
            'R10':0.53134895
            }

np.array([1., 0.24072303, 0.21098072, 0.2211775 , 0.30620438, 0.53134895])

folderstart = "C:/Users/jansenc3/surfdrive/"
folderstart = "C:/Users/Werk/Surfdrive/"
folder = folderstart+"DATA/Power and wavelength data/"
simulationsfolder = folderstart+'DATA/Laser/Simulations/CO2/220517/'
simulationsfolder = folderstart+'DATA/Laser/Simulations/CO2/220603/'

months = {'01':'Jan', '02':'Feb', '03':'Mar', '04':'Apr', '05':'May', '06':'Jun', '07':'Jul', '08':'Aug', '09':'Sep', '10':'Oct', '11':'Nov', '12':'Dec'}
transitions = ['R0', 'R2', 'R4', 'R6', 'R8', 'R10']
# transitions = ['R2']
save=True

def main():
    print('main')

    for transition in transitions:
        datadic = read(transition)
        correct_offsets(transition, datadic)
        # datadic['power'] = absolute_laser_power(datadic['power']) #corrects offset twice, change that
        # plt.plot(datadic['power'], datadic['lockin'])
        # plt.plot(datadic['simulation_power'], datadic['simulation'])
        # plt.show()
        # plt.close()

        # A = fit_amplitude(datadic)
    
        # plt.plot(datadic['power'], A*datadic['lockin'], '.')
        # plt.plot(datadic['simulation_power'], datadic['simulation'])
        # plt.show()
        # plt.close()

        x_guess = np.array([datainfo[transition]['min_power'],1.2,A_guess[transition]])
        pmin, w_power, A = fit_full(datadic, x_guess, bounds=np.array([2,1.5,5]))
        print (pmin, w_power, A)

        plt.plot((datadic['power_raw']-pmin)*w_power, A*(datadic['lockin_raw']-datadic['offresonance']), '.')
        plt.plot(datadic['simulation_power'], datadic['simulation'])
        plt.title(transition+', pmin='+str(np.round(pmin,2))+', w_power='+ str(np.round(w_power,2))+', A='+str(np.round(A,2)))
        plt.xlabel('Power (mW)')
        plt.ylabel('Excited population')
        if save:
            plt.savefig(simulationsfolder+'figures/'+transition+'_fit.png')
        plt.show()
        plt.close()




def fit_full(datadic, x_guess, bounds=None):
    """
    Fits pmin, w_power (measure for laser extinction in window etc), A (measure for how much signal corresponds to population inversion)
    """
    
    bmin = x_guess/bounds
    bmax = x_guess*bounds
    print(bmin,bmax)
    popt = least_squares(fitfunction_full, x_guess, args=([datadic]),bounds=(bmin,bmax))
    return popt.x

def fitfunction_full(x, datadic):
    """
    x = [pmin, w_power, A_experiment]
    """
    pmin = x[0]
    w_power = [1]
    A_experiment = x[2]
    power = (datadic['power_raw'] - pmin)*w_power
    lockin = datadic['lockin_raw'] - datadic['offresonance']

    index = (power > np.min(datadic['simulation_power'])) * (power < np.max(datadic['simulation_power']))
    power_for_fit = power[index]
    data_for_fit = lockin[index]
    sim_interp = np.interp(power_for_fit, datadic['simulation_power'], datadic['simulation'])

    return A_experiment*data_for_fit - sim_interp


def fit_amplitude(datadic, A_experiment_guess=1):
    index = (datadic['power'] > np.min(datadic['simulation_power'])) * (datadic['power'] < np.max(datadic['simulation_power']))
    power_for_fit = datadic['power'][index]
    data_for_fit = datadic['lockin'][index]
    sim_interp = np.interp(power_for_fit, datadic['simulation_power'], datadic['simulation'])


    popt = least_squares(fitfunction_amplitude, np.array([A_experiment_guess]), args=(data_for_fit,sim_interp))
    return popt.x[0]


def fitfunction_amplitude(x, experiment, simulation):
    """
    Returns differences between experiment and simulation that will be minimized.
    Experiment and simulation should correspond to the same power array
    Experiment is the measured data without the offresonance signal
    """
    A_experiment = x[0]
    return simulation - A_experiment*experiment

def correct_offsets(transition, datadic):
    """
    Returns measured power and signal without offsets and power in mW
    """
    datadic['power'] = datadic['power_raw']-datainfo[transition]['min_power']
    datadic['lockin'] = datadic['lockin_raw']-datadic['offresonance']

def read(transition):
    info = datainfo[transition]
    folder = (folderstart+'DATA/'+info['year']+'/'+info['month']+' '+months[info['month']]+'/'+info['year'][-2:]+info['month']+info['day']+'/Laser/Images/'+transition+'/')    
    
    dic = {}
    dic['power_raw'], dic['lockin_raw'] = np.loadtxt(folder+'lockin_vs_power.txt', unpack=True)
    dic['power_raw'] *= 1000
    dic['offresonance'] = np.loadtxt(folder+'lockin_offresonance.txt')

    dic['simulation_power'], dic['simulation'] = np.loadtxt(simulationsfolder+transition+'_averaged.txt', unpack=True)
    # dic['simulation_power'], dic['simulation'] = np.loadtxt(simulationsfolder+transition+'/averaged.txt', unpack=True)
    return dic


def absolute_laser_power(power, T_air=0.90, T_window=0.85): #power_before in mW, power_loss_air in %
    """
    T_air: transmission of the air pocket between the window and power meter \n
    T_window: transmission of the IR window
    """
    return (power-np.min(power))/T_air/T_window

main()