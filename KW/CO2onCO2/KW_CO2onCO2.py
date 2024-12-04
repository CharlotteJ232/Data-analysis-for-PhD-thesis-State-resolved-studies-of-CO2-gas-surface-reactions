import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
import datetime as dt
import os
# import colorednoise as cn

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

# 85K, with laser modulation measured
folder = '2021/11 Nov/211104/KW/'
file = 'kw01'
onresonance = True
modulation_frequency = 1 #Hz
start = 200
stop = 700

typ = 'lin' #for fitting
measured_laser = True
laserfile = 'lasermodulation01'
max_timeshift = 100

# # 85K, with laser modulation measured
# folder = '2021/11 Nov/211104/KW/'
# file = 'kw02'
# onresonance = True
# modulation_frequency = 2 #Hz
# start = 200
# stop = 800

# typ = 'log' #for fitting
# b_guess = -0.01 #for fitting, change if it does not work
# measured_laser = True
# laserfile = 'lasermodulation02'
# max_timeshift = 1.5

# 85K, with laser modulation measured
folder = '2021/11 Nov/211105/KW/'
file = 'kw01'
onresonance = True
modulation_frequency = 2 #Hz
start = 100 #cannot be 0 for some reason
stop = 1700

measured_laser = True
max_timeshift = 2
UTI = True

# 85K, laser not entering the chamber
folder = '2021/11 Nov/211108/KW/'
file = 'kw02'
onresonance = False
modulation_frequency = 2 #Hz
start = 100 #cannot be 0 for some reason
stop = 1450

measured_laser = True
max_timeshift = 2
UTI = True
# 85K, laser off resonance but entering chamber
file = 'kw03'
# # beam not entering chamber, laser off
# file = 'test01'
# # beam not entering chamber, laser off resonance but entering chamber
# file = 'test02'

# 95K, no laser, testing UTI pc measurement
folder = '2021/11 Nov/211111/KW/'
file = 'KW04'
onresonance = False
modulation_frequency = 3 #Hz
start = 60 #cannot be 0 for some reason
stop = 500

measured_laser = False
max_timeshift = 2
UTI = True

# 85K
folder = '2021/11 Nov/211112/KW/'
onresonance = False
modulation_frequency = 3 #Hz
start = 50 #cannot be 0 for some reason
stop = 850

measured_laser = True
max_timeshift = 2
UTI = True
#Laser on resonance, damping
file = 'KW03'
#laser on resonance, no damping
file = 'KW02'
#laser off, no damping
file = 'KW01'
#laser off resonance, damping
file = 'KW04'
# laser on resonance, damping
file = 'KW05'

# 85K
folder = '2021/11 Nov/211115/KW/'
onresonance = False
modulation_frequency = 5 #Hz
start = 50 #cannot be 0 for some reason
stop = 850
measured_laser = True
max_timeshift = 1
UTI = True
#Laser on resonance, damping
file = 'KW01'
#Laser off resonance, damping
file = 'KW02'

# On cryostat, ~77K
folder = '2022/01 Jan/220114/KW/'
onresonance = False
modulation_frequency = 5 #Hz
start = 75 #cannot be 0 for some reason
stop = 450
measured_laser = True
max_timeshift = 1
UTI = True
#test measurement for noise spectrum of the new uti controller
file = 'test'
#off and on resonance???
file = 'KW01'
#off resonance
file = 'KW02'
# #on resonance
file = 'KW03'

# #Noise spectrum test of the new pressure gauge
# folder = '2022/03 Mar/220307/Pressuregauge/'
# onresonance = False
# modulation_frequency = 5 #Hz
# start = 1 #cannot be 0 for some reason
# stop = 450
# measured_laser = True
# max_timeshift = 1
# UTI = True
# #test measurement for noise spectrum of the new uti controller
# file = 'Without_shielding'

# #Noise spectrum test of the new pressure gauge
# folder = '2022/04 Apr/220404/KW/'
# onresonance = False
# modulation_frequency = 5 #Hz
# start = 1 #cannot be 0 for some reason
# stop = 450
# measured_laser = True
# max_timeshift = 1
# UTI = True
# #test measurement for noise spectrum of the new uti controller
# file = 'test_FFT3'


################ General Parameters #################

timediff = 2082844800 #difference between 1-1-1904 (labview) and 1-1-1970 (everything else)
folderstart = "C:/Users/jansenc3/surfdrive/DATA/"
# folderstart = "C:/Users/Werk/Surfdrive/DATA/"
ext = '.txt'
subsavefolder = '' #include / !
savefolder = folderstart+folder+file+'/'+subsavefolder
plot_all=True
save_all=True
use_dummydata = False

#Script
def main():
    if UTI:
        if measured_laser:
            time, data_laser, data = np.loadtxt(folderstart+folder+file+ext,skiprows=3,unpack=True)
        else:
            time, data = np.loadtxt(folderstart+folder+file+ext,skiprows=3,unpack=True)
        time -= time[0]



        plt.plot(time,data)
        if save_all:
            if not os.path.exists(savefolder):
                os.makedirs(savefolder)
            plt.savefig(savefolder+'all_data.png',dpi=500)
        plt.show()
        plt.close()

        if not measured_laser:
            time, data = remove_UTInoise(time, data, factor=1.2)

        #randomize array
        # rng = np.random.default_rng()
        # rng.shuffle(data)
        # # data += 0.001*time
        # # data = 2+4*np.random.rand(len(time))


        time_cut, data = remove_beginend(time, data, start=start,stop=stop, plot=False) 
        plt.plot(time_cut,data)
        if save_all:
            plt.savefig(savefolder+'cut_data.png',dpi=500)
        plt.show()
        plt.close()


        #CALCULATE FFT
        samplingrate = 100 #Hz
        artificialsampling = np.arange(np.min(time_cut), np.max(time_cut), 1/samplingrate) 
        data_fft = np.interp(artificialsampling, time_cut, data)

        fftfreq = np.fft.rfftfreq(len(data_fft), 1/samplingrate)
        fft = np.fft.rfft(data_fft)

        np.savetxt(savefolder+'fft_'+str(start)+'_'+str(stop)+'.txt', np.column_stack((fftfreq.astype(float), np.absolute(fft))), header='Freq, y')
        plt.plot(fftfreq, np.abs(fft))
        x_min = 0
        x_max =50
        plt.axis([x_min, x_max,0,5])
        plt.grid()
        plt.xlabel('Freq (Hz)')
        plt.title('Signal FFT')
        if save_all:
            plt.savefig(savefolder+'fft'+str(x_min)+'_'+str(x_max)+'.png',dpi=500)
        plt.show()
        plt.close()

        if measured_laser:
            data_laser_fft = np.interp(artificialsampling, time, data_laser)
            fft = np.fft.rfft(data_laser_fft)
            fftfreq = np.fft.rfftfreq(len(data_laser_fft), 1/samplingrate)
            np.savetxt(savefolder+'fft_laser_'+str(start)+'_'+str(stop)+'.txt', np.column_stack((fftfreq.astype(float), np.absolute(fft))), header='Freq, y')
            plt.plot(fftfreq, np.abs(fft))
            plt.grid()
            plt.axis([x_min, x_max,0,100])
            plt.xlabel('Freq (Hz)')
            plt.title('Laser FFT')
            if save_all:
                plt.savefig(savefolder+'fft_laser'+str(x_min)+'_'+str(x_max)+'.png',dpi=500)
            plt.show()
            plt.close()



        #CALCULATE DIFFERENCE
        if measured_laser:
            difference_array = []
            timeshiftarray = np.arange(-max_timeshift,max_timeshift,0.01)
            filenameaddition = ''
            if use_dummydata:
                #input values
                beta = 1         # the exponent: 0=white noite; 1=pink noise;  2=red noise (also "brownian noise")
                samples = len(time_cut)  # number of samples to generate (time series extension)
                #Deffing some colores
                data = cn.powerlaw_psd_gaussian(beta, samples)
                filenameaddition = str(data[0]) #random number so figure is not overwritten when rerunning the code
                data = dummydata(time_cut)

            for timeshift in timeshiftarray:
                # bool_laser has the same length as time_cut, but a different part of the laser modulation signal is selected each time by convert_bool_lasermodulation.
                # It is used as a data mask for the same dataset each time. Just the shape of the data mask is different. 
                bool_laser = convert_bool_lasermodulation(time, data_laser, time_cut+timeshift)
                # print (len(data[bool_laser]), len(data[np.invert(bool_laser)]))
                difference = np.average(data[bool_laser]) - np.average(data[np.invert(bool_laser)])
                difference_array.append(difference)
            np.savetxt(savefolder+'timeshift_'+str(start)+'_'+str(stop)+filenameaddition+'.txt', np.column_stack((timeshiftarray,difference_array)), header='timeshift, difference (V)')
            plt.plot(timeshiftarray, difference_array)
            plt.plot([-timeshift, timeshift],[0,0], color='black')
            plt.grid(linestyle='--')
            plt.xlabel('Timeshift (s)')
            plt.ylabel('difference on/off resonance')
            if save_all:
                plt.savefig(savefolder+'difference_'+str(start)+'-'+str(stop)+'_'+str(max_timeshift)+filenameaddition+'.png',dpi=500)
            plt.show()
            plt.close()



#for older data only
    else:
        time, data = np.loadtxt(folderstart+folder+file+ext,skiprows=3,unpack=True)
        with open(folderstart+folder+file+ext, 'r') as qs_file:
            starttime = float(qs_file.readlines()[1])
        if measured_laser:
            time_laser, data_laser = np.loadtxt(folderstart+folder+laserfile+ext,unpack=True)
            time_laser -= timediff
        plt.plot(time,data)
        if save_all:
            if not os.path.exists(savefolder):
                os.makedirs(savefolder)
            plt.savefig(savefolder+'all_data.png',dpi=500)
        plt.show()
        plt.close()
        time, data = remove_beginend(time, data, start=start,stop=stop, plot=False)
        level_data = remove_background(time,data,typ=typ)

        difference_array = []
        timeshiftarray = np.arange(-max_timeshift,max_timeshift,0.01)
        for timeshift in timeshiftarray:
            bool_laser = convert_bool_lasermodulation(time_laser, data_laser, time+starttime+timeshift)
            difference = np.average(level_data[bool_laser]) - np.average(level_data[np.invert(bool_laser)])
            difference_array.append(difference)
            # plt.plot(time[bool_laser], level_data[bool_laser], '.')
            # plt.plot(time[np.invert(bool_laser)], level_data[np.invert(bool_laser)], '.')
            # plt.show()
            # plt.close()
        plt.plot(timeshiftarray, difference_array)
        plt.xlabel('Timeshift (s)')
        plt.ylabel('difference on/off resonance')
        if save_all:
            plt.savefig(savefolder+'difference_'+str(start)+'-'+str(stop)+'_'+str(max_timeshift)+'.png',dpi=500)
        plt.show()
        plt.close()
    

        #BELOW HERE IS THE OLD CODE, WHICH IS NOT USED ANYMORE
        # stacked_time = stack_data(time, modulation_frequency)
        
        # plt.plot(stacked_time, level_data, '.')
        # plt.plot(stacked_time,np.full(len(stacked_time),np.average(level_data)))
        # if save_all:
        #     plt.savefig(savefolder+'stacked_data.png',dpi=500)
        # plt.show()
        # plt.close()

        # averaged_data, std = average_data(stacked_time,level_data)
        # plt.errorbar(np.unique(stacked_time), averaged_data, yerr=std, capsize=5, marker='o',ls='')
        # if save_all:
        #     plt.savefig(savefolder+'averaged_data.png',dpi=500)
        # plt.show()
        # plt.close()

        # find_wave(stacked_time, level_data)


    print('the end')

def dummydata(time):
    data = np.zeros(len(time))
    for freq in np.arange(0.1,100,0.1):
        phase = np.random.random()*2*np.pi
        phase = 0
        data += 1/freq*np.sin(2*np.pi*freq*time+phase)
    return data

def whitenoise(time):
    return np.random.randn(len(time))


def remove_UTInoise(time, data, factor=2):
    remove1 = np.argwhere(data > factor*np.roll(data,1))
    remove2 = np.argwhere(data > factor*np.roll(data,-1))
    remove = np.concatenate((remove1, remove2))
    return np.delete(time, remove), np.delete(data, remove)


    
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


def convert_bool_lasermodulation(time_laser, data_laser, time_kw):
    interpolated_laser = np.interp(time_kw, time_laser, data_laser)
    thresh = np.average(data_laser)
    bool_laser = interpolated_laser > thresh
    return bool_laser

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
    """
    This is not used anymore!
    """
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








