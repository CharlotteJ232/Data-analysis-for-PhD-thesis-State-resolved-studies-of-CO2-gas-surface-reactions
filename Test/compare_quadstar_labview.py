import numpy as np
from matplotlib import pyplot as plt

folderstart = 'C:/Users/jansenc3/surfdrive/Data/2021/'
# folderstart = 'C:/Users/Werk/surfdrive/DATA/2021/'
folder = folderstart+'11 Nov/211104/Testing/'

labview_filename = 'labview01.txt'
quadstar_filename = 'quadstar03.txt'

labview_path = folder+labview_filename
# labview_path = labview_filename
quadstar_path = folder+quadstar_filename
# quadstar_path = quadstar_filename

timediff = 2082844800 #difference between 1-1-1904 (labview) and 1-1-1970 (everything else)

def main():
    labview_time, labview_data = np.loadtxt(labview_path, unpack=True)
    quadstar_time, quadstar_data = np.loadtxt(quadstar_path, skiprows=3, unpack=True)
    with open(quadstar_path, 'r') as qs_file:
        quadstar_starttime = float(qs_file.readlines()[1])
    quadstar_time += quadstar_starttime
    quadstar_time += timediff
    quadstar_data /= np.max(quadstar_data)
    labview_data /= np.max(labview_data)

    plt.plot(labview_time, labview_data, label='LabView')
    plt.plot(quadstar_time, quadstar_data, label='Quadstar')
    plt.axis([np.min(quadstar_time), np.min(quadstar_time)+10,None,None])
    plt.legend()
    plt.xlabel('timestamp')
    plt.ylabel('normalized signals')
    plt.savefig(folder+'Images/'+labview_filename+quadstar_filename+'.png')
    plt.show()
    plt.close()

main()