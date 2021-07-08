import numpy as np
from matplotlib import pyplot as plt

"""
Currently for stepped (111) only
"""
surface = '111'
cut_angle = 10  #degrees
lattice_parameter = 1
Ebeam_eV = 150 #eV
Ebeam_J = Ebeam_ev * 1.602176634E-19


h = 6.62607015E-34 #Js
m_e = 9.1093837015E-31 #kg
p_e = np.sqrt(2*m_e*Ebeam_J)
wl_electron = h/p_e
Lterrace = hstep / np.tan(cut_angle/180*np.pi)

surfaces = {'111': {'hstep':np.sqrt(1/3)*lattice_parameter,'a':np.sqrt(1/2)*lattice_parameter, 'b':np.sqrt(1/2)*lattice_parameter, 'angle':60*np.pi/180}}

def main():

    plot_lattice(plot=True)

def plot_lattice(n_a=5, n_b=5, plot=False):
    print('not finished')

    
    points_x = []
    points_y = []
    for i_a in range(n_a+1):
        for i_b in range(n_b+1):
            points_x.append(i_a*surfaces[surface]['a']+i_b%2*np.cos(surfaces[surface]['angle'])*surfaces[surface]['b']) #i_b%2 to create square instead of diamond
            points_y.append(i_b*surfaces[surface]['b']) #np.sin(surfaces[surface]['angle'])*
    points_x = np.array(points_x)
    points_y = np.array(points_y)

    plt.plot(points_x,points_y, 'o')
    plt.title("Real space")

    return points_x, points_y        



def plot_reciprocal():
    print('not finished')

main()