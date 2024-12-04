import numpy as np
import glob
from matplotlib import pyplot as plt
from matplotlib import cm 
import datetime as dt
from scipy.optimize import curve_fit, least_squares, minimize
from joblib import Parallel, delayed

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
            } #note: for fitting, the inverse of this is used

A_guess = { 'R0':3.0,
            'R2':1,
            'R4':0.8,
            'R6':0.8,
            'R8':1.2,
            'R10':2
            } #note: for fitting, the inverse of this is used



folderstart = "C:/Users/jansenc3/surfdrive/"
folderstart = "C:/Users/Werk/Surfdrive/"
folder = folderstart+"DATA/Power and wavelength data/"
simulationsfolder = folderstart+'DATA/Laser/Simulations/CO2/220602/' #without rotation
simulationsfolder = folderstart+'DATA/Laser/Simulations/CO2/220603/' #with rotation
simulationsfolder = folderstart+'DATA/Laser/Simulations/CO2/230817/withrotations/' #with rotation
# simulationsfolder = folderstart+'DATA/Laser/Simulations/CO2/231206/test/' #without rotation and projection

months = {'01':'Jan', '02':'Feb', '03':'Mar', '04':'Apr', '05':'May', '06':'Jun', '07':'Jul', '08':'Aug', '09':'Sep', '10':'Oct', '11':'Nov', '12':'Dec'}
transitions = ['R0', 'R2', 'R4', 'R6', 'R8', 'R10']
# transitions = ['R2']
save=False

def main():
    print('main')

#### Fit with separate width ####
    # amplitudes = []
    # for transition in transitions:
    #     datadic = read(transition)

    #     x_guess = np.array([1.2,datainfo[transition]['min_power'],1/A_guess[transition]/3]) #1/A_guess/3 because I changed the fit parameter and I divided the data by the off-resonance signal (close to 3 V)
    #     print(x_guess)
    #     w_power, pmin, A = fit_full(datadic, x_guess, bounds=np.array([1.1,1.5,3]))
    #     if save:
    #         np.savetxt(simulationsfolder+transition+'_x.txt', [w_power, pmin, A])
    #     print (pmin, w_power, A)
    #     amplitudes.append(A)
    #     plot_fit(transition, datadic, w_power, pmin, A, name='')
    # print(amplitudes)

#### Fit with shared width ####
    popt = fit_all()
    # print(popt)
    x = popt.x
    hess = popt.hess
    print(1/np.diag(hess))
    if save:
        np.savetxt(simulationsfolder+'x.txt', x)
    print(popt)
    # x = [1.09090909, 4.86848057, 0.98337631, 4.32269159, 4.15274705, 3.80092913, 5.25014092, 5.10166111, 4.64286135, 4.88022767, 3.1502073, 4.4786131, 1.94651]
    index = 1
    for transition in transitions:
        datadic = read(transition)
        w_power = x[0]
        pmin = x[index]
        A = x[index+1]
        index += 2
        plot_fit(transition, datadic, w_power, pmin, A, name='_all')




def plot_fit(transition, datadic, w_power, pmin, A, name=''):
        plt.plot((datadic['power_raw']-pmin)*w_power, (datadic['lockin_raw']-datadic['offresonance'])/datadic['offresonance'], '.')
        plt.plot(datadic['simulation_power'], datadic['simulation']*A)
        plt.title(transition+', pmin='+str(np.round(pmin,3))+', w_power='+ str(np.round(w_power,3))+', A='+str(np.round(A,3)))
        plt.xlabel('Power (mW)')
        plt.ylabel('Excited population')
        if save:
            if not os.path.exists(simulationsfolder+'figures/'):
                os.makedirs(simulationsfolder+'figures/')
            plt.savefig(simulationsfolder+'figures/'+transition+'_fit'+name+'.png')
        plt.show()
        plt.close()

def fit_all(w_power_guess = 1.3, bounds=[1.2,1.5,3]):
    """
    x_guess: w_power, pmin, A, pmin, A, etc
    """
    x_guess = [w_power_guess]
    bounds_for_minimize = [(w_power_guess/bounds[0], w_power_guess*bounds[0])]
    for transition in transitions:
        pmin = datainfo[transition]['min_power']
        x_guess.append(pmin)
        bounds_for_minimize.append((pmin/bounds[1],pmin*bounds[1]))
        A = 1/A_guess[transition]
        x_guess.append(A)
        bounds_for_minimize.append((A/bounds[2],A*bounds[2]))
    popt = minimize(fitfunction_all, x_guess,bounds=bounds_for_minimize)
    return popt


def fitfunction_all(x):
    sum = 0
    index = 1
    # num_cores = min(cpu_count()-1, len(transitions))
    # Parallel(n_jobs=num_cores)(delayed(simulate_trajectory)() for transition in transitions)
    for transition in transitions:
        datadic = read(transition)
        sum += fitfunction_full([x[0], x[index], x[index+1]], datadic)
        index += 2
    return sum

# def fitfunction_all_loop(transition, x):

def fit_full(datadic, x_guess, bounds=None):
    """
    Fits pmin, w_power (measure for laser extinction in window etc), A (measure for how much signal corresponds to population inversion)
    """
    bmin = x_guess/bounds
    bmax = x_guess*bounds
    bounds = [(bmin[i], bmax[i]) for i in range(len(x_guess))]
    popt = minimize(fitfunction_full, x_guess, args=(datadic),bounds=bounds)
    # popt = least_squares(fitfunction_full, x_guess, args=([datadic]))
    return popt.x

def fitfunction_full(x, datadic):
    """
    x = [pmin, w_power, A_experiment]
    """
    w_power = x[0]
    pmin = x[1]
    A_experiment = x[2]
    power = (datadic['power_raw'] - pmin)*w_power
    lockin = (datadic['lockin_raw'] - datadic['offresonance'])/datadic['offresonance'] #division by offresonance to correct for changes in sensitivity

    index = (power > np.min(datadic['simulation_power'])) * (power < np.max(datadic['simulation_power']))
    power_for_fit = power[index]
    data_for_fit = lockin[index]
    sim_interp = np.interp(power_for_fit, datadic['simulation_power'], datadic['simulation'])

    # print(np.shape(power_for_fit), np.shape(data_for_fit), np.shape(sim_interp))
    return np.sum((data_for_fit/A_experiment - sim_interp)**2)/len(sim_interp)


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

    # dic['simulation_power'], dic['simulation'] = np.loadtxt(simulationsfolder+transition+'_averaged.txt', unpack=True)

    filename = glob.glob(simulationsfolder+transition+'*.txt')
    # print(filename)
    dic['simulation_power'], dic['simulation'] = np.loadtxt(filename[0], unpack=True)

    # dic['simulation_power'], dic['simulation'] = np.loadtxt(simulationsfolder+transition+'/averaged.txt', unpack=True)
    return dic


def absolute_laser_power(power, T_air=0.90, T_window=0.85): #power_before in mW, power_loss_air in %
    """
    T_air: transmission of the air pocket between the window and power meter \n
    T_window: transmission of the IR window
    """
    return (power-np.min(power))/T_air/T_window

main()