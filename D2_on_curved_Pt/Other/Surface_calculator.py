# -*- coding: utf-8 -*-
"""
Created on Fri Feb 28 15:05:09 2020
NOTE: this calculates probed step density with convolution of the molecular beam. 
@author: jansenc3
"""

import numpy as np
from matplotlib import pyplot as plt
plt.style.reload_library()
#plt.style.use('voorbeeld')

step_size = 0.001 #mm, step size for the positions array. nothing physical, so can be chosen
positions_array = np.arange(-14,14,step_size) #positions in mm
crystal_radius = 15 #mm, radius of curvature of the crystal
step_height = 0.2265637926 #nm, for steps on pt111
atom_radius = 0.13874 #nm, for platinum atom
row_width = np.sqrt(3)*atom_radius
beam_width = 0.124 #mm


# position = np.absolute(19.9825-22) #mm, position that you want to use for calculations
# angle = np.arcsin(position/crystal_radius)
angle = 5.6 #degrees
angle_rad = angle * np.pi / 180
position = crystal_radius * np.sin(angle_rad)


def main():    
    print ('step density at angle '+str(angle)+' degrees and position '+str(position)+' mm calculated from angle is: '+str(calc_step_density_from_angle(angle_rad, step_height)))
    x_pos=np.arange(-10, 10, 0.001)
    plt.plot(x_pos, np.arcsin(x_pos/crystal_radius))
    plt.show()
    plt.close()
    print ('step density at angle '+str(angle)+' degrees and position '+str(position)+' mm is: '+str(calc_step_density(position, step_height)))

    positions, step_densities = probed_step_density(positions_array)
    
    step_density = np.interp(position, positions, step_densities)
    print ('probed step density at angle '+str(angle)+' degrees and position '+str(position)+' mm is: '+str(step_density))
    
    terrace_fraction(positions, step_densities, plot=True)
    
    unit_cell_length = 1/np.absolute(step_density)/row_width
    print("unit cell length ",unit_cell_length,' atom rows')
    print ('step size ',calc_step_size(unit_cell_length, 0.4,0.05))

        
def calc_step_size(unit_cell_length, s0, s0_center): #in number of atom rows
    return (s0-s0_center)/(1-s0_center)*unit_cell_length

def terrace_fraction(positions, step_density, plot=False):
    """
    assumes step size of single row
    """
    terrace_ratio = 1/np.absolute(step_density)/np.sqrt(3)/atom_radius - 1
    terrace_fraction = terrace_ratio / (terrace_ratio+1)
    
    if plot:
        #        plt.scatter(self.pos, self.terrace_ratio)
        plt.scatter(positions, terrace_fraction*100)
        plt.title('Percentage of terraces')
        plt.xlabel('position on crystal (mm)')
        plt.ylabel('% terrace')
        plt.show()
        plt.close()
    
def calc_step_density(x, step_height):
    """
    returns step density in /nm, because step_height is in nm.
    """
    return (1/step_height * x / np.sqrt(crystal_radius**2-x**2)) #np.absolute(x)

def calc_step_density_from_angle(theta, step_height):
    return (np.tan(theta)/step_height)    

def probed_step_density(positions, plot=True):
    
    step_density = calc_step_density(positions, step_height)
    
    som = np.cumsum(step_density)
    
    n_elements = np.int(beam_width/step_size/2) #number of array elements for calculating convolution of the beam on the surface
    
    probed_step_density = (np.roll(som,-n_elements)-np.roll(som,n_elements))/(n_elements*2)
    probed_step_density = probed_step_density[n_elements:-n_elements]
    positions = positions[n_elements:-n_elements]
    
    if plot:
        plt.plot(positions, np.absolute(probed_step_density))
        plt.title('Step density')
        plt.xlabel('Position relative to center (mm)')
        plt.ylabel('Step density')
    #        plt.axis([-0.2,0.2,0,0.1])
        plt.show()
        plt.close()
        

    return positions, probed_step_density

main()