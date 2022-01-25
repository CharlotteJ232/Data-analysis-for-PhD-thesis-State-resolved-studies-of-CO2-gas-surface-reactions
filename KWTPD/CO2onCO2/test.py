import numpy as np
import matplotlib.pyplot as plt
import os

# 95K, no laser, testing UTI pc measurement
folder = '2021/11 Nov/211111/KW/'
file = 'kw03'
onresonance = False
modulation_frequency = 2 #Hz
start = 1 #cannot be 0 for some reason
stop = 500
subsavefolder = ''
measured_laser = False
max_timeshift = 2
UTI = True


################ General Parameters #################

timediff = 2082844800 #difference between 1-1-1904 (labview) and 1-1-1970 (everything else)
folderstart = "C:/Users/jansenc3/surfdrive/DATA/"
# folderstart = "C:/Users/Werk/Surfdrive/DATA/"
ext = '.txt'
subsavefolder = ''
savefolder = folderstart+folder+file+'/'+subsavefolder
plot_all=True
save_all=False

#Script
def main():
    if UTI:
        time, data = np.loadtxt(folderstart+folder+file+ext,skiprows=3,unpack=True)
        time -= time[0]

        plt.plot(time,data, '.')
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

        samplingrate = 200 #Hz
        artificialsampling = np.arange(np.min(time_cut), np.max(time_cut), 1/samplingrate) 
        data_fft = np.interp(artificialsampling, time_cut, data)
        data_laser_fft = np.interp(artificialsampling, time, data_laser)

        fft = np.fft.rfft(data_fft)
        fftfreq = np.fft.rfftfreq(len(data_fft), 1/samplingrate)
        plt.plot(fftfreq, np.abs(fft))
        plt.axis([0,5,0,5000])
        plt.xlabel('Freq (Hz)')
        plt.title('Signal FFT')
        if save_all:
            plt.savefig(savefolder+'fft.png',dpi=500)
        plt.show()
        plt.close()


main()