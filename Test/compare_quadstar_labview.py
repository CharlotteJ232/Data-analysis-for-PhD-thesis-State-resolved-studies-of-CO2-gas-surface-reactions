import numpy as np
from matplotlib import pyplot as plt



labview_filename = 'Chopper01.txt'
quadstar_filename = 'kw13.txt'

labview_path = labview_filename
quadstar_path = quadstar_filename

timediff = 66*31556926 #because labview starts in 1904 and everything else in 1970
timediff += 86400 #one day (?)

def main():
    labview_time, labview_data = np.loadtxt(labview_path, unpack=True)
    quadstar_time, quadstar_data = np.loadtxt(quadstar_path, skiprows=3, unpack=True)
    with open(quadstar_path, 'r') as qs_file:
        quadstar_starttime = float(qs_file.readlines()[1])
    quadstar_time += quadstar_starttime
    quadstar_time += timediff

    plt.plot(labview_time, labview_data, label='LabView')
    plt.plot(quadstar_time, quadstar_data, label='Quadstar')
    plt.legend()
    plt.show()
    plt.close()

main()