import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
import datetime as dt

#Measurement parameters
#100K
folder = '2021/10 Oct/211015/KW/'
file = 'kw01'
onresonance = True
modulation_frequency = 1 #Hz
start = 200
stop = 640
b_guess = -0.01 #for exponential fit, change if fit does not work

# 85K
folder = '2021/10 Oct/211018/KW/'
file = 'kw04'
onresonance = True
modulation_frequency = 1 #Hz
start = 200
stop = 1300
subsavefolder = ''
b_guess = -0.001 #for exponential fit, change if fit does not work

start = 200
stop = 500
subsavefolder = 'Firstpart/'

start = 1000
stop = 1300
subsavefolder = 'Lastpart/'

file='kw01'
start=1
stop=100

################ General Parameters #################

folderstart = "C:/Users/jansenc3/surfdrive/DATA/"
# folderstart = "C:/Users/Werk/Surfdrive/DATA/"
ext = '.txt'
savefolder = folderstart+folder+file+'/'+subsavefolder
plot_all=True
save_all=False

#Script
def main():
    time, data = np.loadtxt(folderstart+folder+file+ext,skiprows=3,unpack=True)
    plt.plot(time,data)
    if save_all:
        if not os.path.exists(savefolder):
            os.makedirs(savefolder)
        plt.savefig(savefolder+'all_data.png',dpi=500)
    plt.show()
    plt.close()
    time, data = remove_beginend(time, data, start=start,stop=stop, plot=plot_all)

    level_data = remove_background(time,data,typ='log')
    stacked_time = stack_data(time, modulation_frequency)
    
    plt.plot(stacked_time, level_data, '.')
    plt.plot(stacked_time,np.full(len(stacked_time),np.average(level_data)))
    if save_all:
        plt.savefig(savefolder+'stacked_data.png',dpi=500)
    plt.show()
    plt.close()

    averaged_data, std = average_data(stacked_time,level_data)
    plt.errorbar(np.unique(stacked_time), averaged_data, yerr=std, capsize=5, marker='o',ls='')
    if save_all:
        plt.savefig(savefolder+'averaged_data.png',dpi=500)
    plt.show()
    plt.close()

    find_wave(stacked_time, level_data)


    print('the end')

    
def remove_beginend(time, data, start=None, stop=None, begin_s=10, end_s=2, plot=True):
    """
    start = manual start
    stop = manual stop
    begin_s = extra seconds removed when automatically determining start
    end_s = extra seconds removed when automatically determining end
    """
    if start and stop:
        end = time>start
        begin = time<stop
        total = end*begin

    if plot:
        plt.plot(time[total],data[total])
        if save_all:
            plt.savefig(savefolder+'cut_data.png',dpi=500)
        plt.show()
        plt.close()
    return time[total], data[total]



def remove_background(time, data, typ='log', plot=True):
    print('remove_background')
    print('remove_tilt')
    if typ=='log':
        def exp(time, A,b,c,d):
            return A*(1-np.exp(b*(time-c)))+d
        A_guess = np.max(data)-np.min(data)
        # b_guess = -0.01
        c_guess = time[0]
        d_guess = np.min(data)
        guess = [A_guess,b_guess,c_guess,d_guess]
        popt = curve_fit(exp, time,data,p0=guess)
        if plot:
            plt.plot(time, data)
            plt.plot(time, exp(time,*popt[0]))
            if save_all:
                plt.savefig(savefolder+'fitted_data.png',dpi=500)
            plt.show()
            plt.close()
        print(popt[0])
        return data - exp(time,*popt[0])
        
    elif typ=='lin':    
        a, b = np.polyfit(time, data, 1) #ax + b
        if plot:
            plt.plot(time, data)
            plt.plot(time, a*time+b)
            if save_all:
                plt.savefig(savefolder+'fitted_data.png',dpi=500)
            plt.show()
            plt.close()
        return data - (a*time+b)


def stack_data(time, freq):
    print('stack_data')
    return time % (1/freq)

def average_data(stacked_time, level_data):
    times = np.unique(stacked_time)
    averaged_data = []
    std = []
    for time in times:
        args = np.argwhere(stacked_time==time)
        averaged_data.append(np.average(level_data[args]))
        std.append(np.std(level_data[args]))
    averaged_data = np.array(averaged_data)
    return averaged_data, np.array(std)

def find_wave_wrong(stacked_time, averaged_data):
    difference = []
    half_length = int(len(averaged_data)/2)
    for i in range(len(averaged_data)):
        shifted_data = np.roll(averaged_data,i)
        diff = np.average(shifted_data[half_length:]) - np.average(shifted_data[:half_length])
        difference.append(diff)
    plt.plot(np.unique(stacked_time), difference, '.')
    plt.ylabel('Difference between first and second half')
    plt.xlabel('Time shift')
    plt.show()
    plt.close()

def find_wave(stacked_time, level_data, timestep=0.005):
    args = np.argsort(stacked_time)
    sorted_time = stacked_time[args]
    #make the array twice as long to remove boundary
    sorted_time_2 = np.concatenate((sorted_time,sorted_time+1/modulation_frequency))
    sorted_data = level_data[args]
    sorted_data_2 = np.concatenate((sorted_data,sorted_data))

    times = np.arange(0,1/modulation_frequency,timestep)

    difference = []
    T = 1/modulation_frequency
    for timeshift in times:
        indices_low_1 = sorted_time_2 > timeshift
        indices_low_2 = sorted_time_2 < timeshift+T/2
        indices_low = indices_low_1*indices_low_2
        indices_high_1 = sorted_time_2 > timeshift+T/2
        indices_high_2 = sorted_time_2 < timeshift+T
        indices_high = indices_high_1*indices_high_2
        diff = np.average(sorted_data_2[indices_low]) - np.average(sorted_data_2[indices_high])
        difference.append(diff)
    plt.plot(times, difference, '.')
    plt.ylabel('Difference between first and second half of period')
    plt.xlabel('Time shift(s)')
    plt.axis([0,T/2,None,None])
    if save_all:
        plt.savefig(savefolder+'difference.png',dpi=500)
    plt.show()
    plt.close()







main()








