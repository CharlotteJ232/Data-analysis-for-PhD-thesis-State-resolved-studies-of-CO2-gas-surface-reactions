from dataclasses import dataclass
from matplotlib import pyplot as plt
import numpy as np
from scipy.optimize import least_squares
from matplotlib.ticker import (AutoMinorLocator, MultipleLocator)
import os
from matplotlib import cm
import colorcet as cc

folderstart = "C:/Users/jansenc3/surfdrive/DATA/"
# folderstart = "C:/Users/Werk/Surfdrive/DATA/"

datalist = []

# # 20221208 With laser off resonance
# datalist.append({'folder': '2022/12 Dec/221208/KW/',
# 'filename': 'tpd01.txt',
# 'gases':['CO2', 'CO'],
# 'tstart':45,
# 'k/s':1,
# 'temps':[100,400]
# })

# # 20221212 With laser on resonance
# datalist.append({'folder': '2022/12 Dec/221212/KW/',
# 'filename': 'tpd01.txt',
# 'gases':['CO2', 'CO', 'H2O', 'O2'],
# 'tstart':17,
# 'k/s':1,
# 'temps':[100,400]
# })

# # 20221213 With laser off resonance
# datalist.append({'folder': '2022/12 Dec/221213/KW/',
# 'filename': 'tpd01.txt',
# 'gases':['CO2', 'CO', 'H2O', 'O2'],
# 'tstart':28,
# 'k/s':1,
# 'temps':[100,400]
# })

# # 20221216 With laser off resonance
# datalist.append({'folder': '2022/12 Dec/221216/KW/',
# 'filename': 'tpd01.txt',
# 'gases':['CO2', 'CO', 'H2O', 'O2'],
# 'tstart':28,
# 'k/s':1,
# 'temps':[100,400]
# })

# # 20221219 With laser on resonance
# datalist.append({'folder': '2022/12 Dec/221219/KW/',
# 'filename': 'tpd01.txt',
# 'gases':['CO2', 'CO', 'H2O', 'O2'],
# 'tstart':28,
# 'k/s':1,
# 'temps':[100,400]
# })

# 20230127 TPD of CO with a small precoverage of O
datalist.append({'folder': '2023/01 Jan/230127/KW/',
'filename': 'tpd01.txt',
'gases':['CO2', 'CO', 'H2O', 'O2'],
'tstart':120,
'k/s':1,
'temps':[80,300]
})


def main():
    for dataset in datalist:
        dataset['data']={}
        all_data = np.loadtxt(folderstart+dataset['folder']+dataset['filename'], skiprows=3, unpack=False)
        # print(np.shape(all_data))
        dataset['data']['time']=all_data[:,0]
        dataset['data']['temp']=(dataset['data']['time']-dataset['tstart'])*dataset['k/s']+dataset['temps'][0]
        for gas, index in zip(dataset['gases'], range(len(all_data)-1)):
            dataset['data'][gas]=all_data[:,index+1]

        for gas in dataset['gases']:
            plt.plot(dataset['data']['temp'], dataset['data'][gas], label=gas)
            plt.ylabel('QMS current (A)')
            plt.xlabel('Time (s)')
        plt.title('TPD')
        plt.axis([None, None, 1E-9, 3E-9])
        plt.axis([dataset['temps'][0], dataset['temps'][1], None, None])
        # plt.title(gas)
        plt.xlabel('Temperature')
        plt.ylabel('QMS current (A)')
        plt.legend(loc='upper right')
        plt.show()
        plt.close()

        datamask = (dataset['data']['time'] > dataset['tstart']+50) * (dataset['data']['time'] < dataset['tstart']+150)
        datamask = np.invert(datamask) * (dataset['data']['time'] > dataset['tstart'])
        dataset['CO2_amp'], dataset['CO2_shift']=fit_CO2_to_CO(dataset, datamask=datamask)

        plt.plot(dataset['data']['temp'], dataset['data']['CO'], label='CO')
        plt.plot(dataset['data']['temp'], dataset['CO2_amp']*dataset['data']['CO2']+dataset['CO2_shift'], label='CO2')
        plt.axis([dataset['temps'][0], dataset['temps'][1], None, None])
        plt.xlabel('Temperature')
        plt.ylabel('QMS current (A)')  
        plt.legend()
        plt.show()
        plt.close()

        plt.plot(dataset['data']['temp'], dataset['data']['CO']-dataset['CO2_amp']*dataset['data']['CO2']-dataset['CO2_shift'])
        plt.axis([dataset['temps'][0], dataset['temps'][1], None, None])
        plt.xlabel('Temperature')
        plt.ylabel('QMS current (A)')
        plt.show()
        plt.close()



def fit_CO2_to_CO(dataset, datamask=None):
    amplitude = [0.5, 2] #must be a list for the least_squares function
    CO2 = dataset['data']['CO2'][datamask]
    CO = dataset['data']['CO'][datamask]

    def CO_minus_CO2(amplitude):
        # print(np.sum(CO-amplitude[0]*CO2))
        return (CO-amplitude[0]*CO2-amplitude[1])*1E10 #multiply by 1E10 because otherwise the function thinks it is optimized immediately, due to small result


    amp_opt = least_squares(CO_minus_CO2, x0=amplitude)

    # print(amp_opt['x'])

    return amp_opt['x'][0], amp_opt['x'][1]



main()
