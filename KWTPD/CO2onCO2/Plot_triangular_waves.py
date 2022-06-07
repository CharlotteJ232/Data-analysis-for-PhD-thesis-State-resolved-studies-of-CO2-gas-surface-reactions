import numpy as np
from matplotlib import pyplot as plt
import scipy.signal

folderstart = "C:/Users/jansenc3/surfdrive/DATA/"
folderstart = "C:/Users/Werk/Surfdrive/DATA/"
folder = '2022/01 Jan/220114/KW/'
file = 'KW03.txt'

def main():

    x = np.arange(0,200,0.01)
    y = np.zeros(len(x))
    # for i in [1,3,5,7,9]:
    #     triangle = scipy.signal.sawtooth(x*i, width=0.5)
    #     plt.plot(x, triangle)
    #     y += triangle
    # plt.show()
    # plt.close()

    # plt.plot(x, y)
    # plt.show()
    # plt.close()


    #test data mask
    data_5hz = np.sin(5*x*2*np.pi)
    datamask_5hz = data_5hz > 0 

    time, modulation = np.loadtxt(folderstart+folder+file, usecols=(0,1), unpack=True)
    for i in np.arange(0,27,0.05):
        phase=np.random.random()*2*np.pi
        difference = apply_datamask(datamask_5hz, np.sin(i*time*2*np.pi+phase))
        # print (i, difference)
        difference_sum += 



    difference_sum = np.zeros(len(x))
    for i in np.arange(0,27,0.05):
        phase=np.random.random()*2*np.pi
        difference = apply_datamask(datamask_5hz, np.sin(i*x*2*np.pi+phase))
        # print (i, difference)
        difference_sum += difference*scipy.signal.sawtooth(i*x*2*np.pi+phase, width=0.5)

    plt.plot(x,difference_sum)
    plt.axis([0,1,None,None])
    plt.show()
    plt.close()

    for i in range(25): 
        print (i, apply_datamask(datamask_5hz, scipy.signal.square(i*x*2*np.pi)))

def apply_datamask(mask, data, round=5):
    return (np.round(np.average(data[mask])-np.average(data[np.invert(mask)]),round))

def convert_bool_lasermodulation(time_laser, data_laser, time_kw):
    interpolated_laser = np.interp(time_kw, time_laser, data_laser)
    thresh = np.average(data_laser)
    bool_laser = interpolated_laser > thresh
    return bool_laser

main()

