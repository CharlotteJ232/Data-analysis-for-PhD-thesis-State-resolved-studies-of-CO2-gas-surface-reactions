import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
import datetime as dt
import os
# import colorednoise as cn

#Test with cold cathode pressure gauge
folder = '2023/05 May/230504/Pressure/'
onresonance = False
modulation_frequency = 5 #Hz
start = 1 #cannot be 0 for some reason
stop = 2100
max_timeshift = 1
#test measurement for noise spectrum of the new uti controller
file = 'Pressure_fast'
file = 'Test2'
measured_laser = False
file_laser = None

#Test with cold cathode pressure gauge
folder = '2023/05 May/230508/Pressure/'
onresonance = False
modulation_frequency = 5 #Hz
start = 750 #cannot be 0 for some reason
stop = 9999999
max_timeshift = 1
#test measurement for noise spectrum of the new uti controller
file = 'background'
measured_laser = False
file_laser = None

# #Test with QMA200
# folder = '2023/05 May/230508/KW/'
# onresonance = False
# modulation_frequency = 3 #Hz
# start = 900 #cannot be 0 for some reason
# stop = 9999999
# max_timeshift = 1
# #test measurement for noise spectrum of the new uti controller
# file = 'KW03'
# measured_laser = False
# file_laser = None

#measurement cold cathode pressure gauge
folder = '2023/05 May/230508/Pressure/'
onresonance = True
modulation_frequency = 3 #Hz
start = 3000 #cannot be 0 for some reason
start = 1
start = 2000
stop = 5200
max_timeshift = 1
file = 'onresonance_R4'
measured_laser = True
file_laser = 'onresonance_R4_laser'

# # measurement cold cathode pressure gauge, 5 Hz
# folder = '2023/05 May/230511/Pressure/'
# onresonance = True
# modulation_frequency = 5 #Hz
# start = 2000 #cannot be 0 for some reason
# start = 1
# stop = 99999
# max_timeshift = 1
# file = 'onresonance'
# measured_laser = False
# file_laser = None

axis=None
axis=[0, 6, 0, 5E-9]
axis_laser=[0, None, 0, 50000]

################ General Parameters #################

timediff = 2082844800 #difference between 1-1-1904 (labview) and 1-1-1970 (everything else)
folderstart = "C:/Users/jansenc3/surfdrive/DATA/"
# folderstart = "C:/Users/Werk/Surfdrive/DATA/"
ext = '.txt'
subsavefolder = '' #include / !
savefolder = folderstart+folder+file+'/'+subsavefolder
plot_all=True
save_all=False
use_dummydata = False

################ Script ############################
def main():
    if measured_laser:
        time_laser, data_laser= np.loadtxt(folderstart+folder+file_laser+ext,skiprows=3,unpack=True, usecols=(0,1))
        time_laser -= time_laser[0]
    time, data = np.loadtxt(folderstart+folder+file+ext,skiprows=3,unpack=True, usecols=(0,1))
    time -= time[0]

    # #for testing the effect of laser going off-resonance
    if use_dummydata:
        data = dummydata(time, 5, peak_datamask=True)

    plt.plot(time,data)
    if save_all:
        if not os.path.exists(savefolder):
            os.makedirs(savefolder)
        plt.savefig(savefolder+'all_data.png',dpi=500)
    plt.show()
    plt.close()

    time_cut, data = remove_beginend(time, data, start=start,stop=stop, plot=False) 
    plt.plot(time_cut,data)
    if save_all:
        plt.savefig(savefolder+'cut_data.png',dpi=500)
    plt.show()
    plt.close()


    #CALCULATE FFT
    samplingrate = 100 #Hz
    samplingrate_from_data = 1/np.average(np.diff(time_cut))
    samplingrate = samplingrate_from_data
    artificialsampling = np.arange(np.min(time_cut), np.max(time_cut), 1/samplingrate) 
    data_fft = np.interp(artificialsampling, time_cut, data)

    fftfreq = np.fft.rfftfreq(len(data_fft), 1/samplingrate)
    fft = np.fft.rfft(data_fft)

    np.savetxt(savefolder+'fft_'+str(start)+'_'+str(stop)+'.txt', np.column_stack((fftfreq.astype(float), np.absolute(fft))), header='Freq, y')
    plt.plot(fftfreq, np.abs(fft))
    x_min = 0
    x_max =samplingrate/2
    plt.axis([x_min, x_max,0,np.max(np.abs(fft[100:]))]) #50: because otherwise the 1/f noise would dominate
    if axis:
        plt.axis(axis)
    plt.grid()
    plt.xlabel('Freq (Hz)')
    plt.title('Signal FFT')
    if save_all:
        plt.savefig(savefolder+'fft'+str(x_min)+'_'+str(x_max)+'.png',dpi=500)
    plt.show()
    plt.close()

    if measured_laser:
        data_laser_fft = np.interp(artificialsampling, time_laser, data_laser)
        fft = np.fft.rfft(data_laser_fft)
        fftfreq = np.fft.rfftfreq(len(data_laser_fft), 1/samplingrate)
        np.savetxt(savefolder+'fft_laser_'+str(start)+'_'+str(stop)+'.txt', np.column_stack((fftfreq.astype(float), np.absolute(fft))), header='Freq, y')
        plt.plot(fftfreq, np.abs(fft))
        plt.grid()
        plt.axis(axis_laser)
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
            beta = 1         # the exponent: 0=white noite; 1=pink / 1/f noise;  2=red noise (also "brownian noise")
            samples = len(time_cut)  # number of samples to generate (time series extension)
            #Deffing some colores
            data = cn.powerlaw_psd_gaussian(beta, samples)
            filenameaddition = str(data[0]) #random number so figure is not overwritten when rerunning the code
            data = dummydata(time_cut)

        for timeshift in timeshiftarray:
            # bool_laser has the same length as time_cut, but a different part of the laser modulation signal is selected each time by convert_bool_lasermodulation.
            # It is used as a data mask for the same dataset each time. Just the shape of the data mask is different. 
            bool_laser = convert_bool_lasermodulation(time_laser, data_laser, time_cut+timeshift)
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


    print('the end')

def dummydata(time, peakf=None, peak_datamask=False):
    """
    Creates 1/f noise with a peak at a certain frequency
    can also use a datamask to simulate laser going off resonance
    Turns out periods of off-resonance makes the peak smaller, a little more
    than proportional to the time spent off-resonance
    """
    data = np.zeros(len(time))
    for freq in np.arange(0.1,100,0.1):
        phase = np.random.random()*2*np.pi
        phase = 0
        data += 1/freq*np.sin(2*np.pi*freq*time+phase)
    if peakf:
        sine = 1/peakf*np.sin(2*np.pi*peakf*time) #at phase 0
        if peak_datamask:
            max = len(data)
            datamask = np.full(max, True)
            datamask[int(max*0.2):int(max*0.3)]=False
            datamask[int(max*0.5):int(max*0.51)]=False
            datamask[int(max*0.6):int(max*0.75)]=False
            datamask[int(max*0.9):int(max*0.95)]=False
            total = 0.1+0.01+0.15+0.05
            sine *= datamask/(1-total)
            plt.plot(time, sine)
            plt.show()
            plt.close()
        
        data += sine
    return data

def whitenoise(time):
    return np.random.randn(len(time))

    
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








