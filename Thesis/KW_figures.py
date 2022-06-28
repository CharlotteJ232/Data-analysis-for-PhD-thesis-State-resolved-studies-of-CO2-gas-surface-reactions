import pickle
import numpy as np
from matplotlib import pyplot as plt
import sys



folderstart = 'P:/Surfdrive/'
folderstart = 'C:/Users/Werk/surfdrive/'

sys.path.append(folderstart+'Python/')

from KWTPD.process_KW_OO_positions import Measurement

folder = folderstart+'DATA/2020/03 Mar/200311/KW/Processed_data/'
positions = np.arange(179, 220, 2)

with open(folder+'data.pickle', 'rb') as f:
    # The protocol version used is detected automatically, so we do not
    # have to specify it.
    measurements = pickle.load(f)


for position in positions:
    measurements[position].plot_measurement_set(save=False)
