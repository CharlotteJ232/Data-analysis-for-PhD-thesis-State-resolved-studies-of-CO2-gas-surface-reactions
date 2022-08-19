from dataclasses import dataclass
from matplotlib import pyplot as plt
import numpy as np
from scipy.optimize import least_squares
from brukeropusreader import read_file
from matplotlib.ticker import (AutoMinorLocator, MultipleLocator)
import os
from matplotlib import cm
import colorcet as cc

folderstart = "C:/Users/jansenc3/surfdrive/DATA/"
folderstart = "C:/Users/Werk/Surfdrive/DATA/"
folder = '2022/07 Jul/220728/RAIRS/'
filestart = '5CO2_15He.00'

@dataclass
class Experiment:
    def add_spectrum():
        print('add_spectrum')

@dataclass
class Spectrum:
    beam_molecules: dict = None #example: {'CO2':5, 'He':15}
    dosing_time: float = None #total dosing time of the molecule of interest
    waiting_time: float = None #total time since the start of the experiment, can be used to account for desorption losses etc.

    def load_data(self, filename):
        self.all_data = read_file(filename)
        #calculate frequency array for plotting
        self.nmin = self.all_data['Fourier Transformation (Rf)']['HFQ']
        self.nmax = self.all_data['Fourier Transformation (Rf)']['LFQ']
        n_datapoints = len(self.all_data['AB'])
        self.wavenumbers = np.linspace(self.nmax, self.nmin, n_datapoints)
        self.absorbance = self.all_data['AB']
        self.background = self.all_data['ScRf']
        self.raw = self.all_data['ScSm']

    def set_datamask(self, datamask):
        self.datamask = datamask

    def calculate_ab_background(self, datamask=None, typ='lin'):
        print('ab_background')
        if not self.datamask:
            if datamask:
                self.set_datamask(datamask)
            else:
                print('Please provide datamask')
        if typ == 'lin':
            a, b = np.polyfit(self.wavenumbers[datamask], self.absorbance[self.datamask], 1)
            self.ab_background=a*self.wavenumbers + b



    def integral(self, nmin, nmax, name=None):
        print('integral')




def main():

    testspectrum = Spectrum()
    testspectrum.load_data(folderstart+folder+filestart+'20')
    # print(testspectrum.all_data.keys())

    plt.plot(testspectrum.wavenumbers, testspectrum.absorbance)
    plt.axis([testspectrum.nmax, testspectrum.nmin, None, None])
    plt.show()
    plt.close()







main()


