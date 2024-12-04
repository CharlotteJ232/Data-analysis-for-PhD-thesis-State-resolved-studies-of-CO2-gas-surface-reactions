import numpy as np
from matplotlib import pyplot as plt

# most calculations based on valleau & deckers 1965

kB = 1.38064852E-23 #boltzmann constant, J/K
R = 8.314 #gas constant per mole

T_nozzle = 300 #Kelvin, temperature in nozzle
d_nozzle = 28E-6 #m, nozzle diameter

pressure_0 = 4970E2 #pa, pressure in nozzle
n_0 = pressure_0 / kB / T_nozzle #/m3, number density in nozzle

d_skim = 0.45E-3 #m, skimmer diameter, a in text
s_skim = 10E-3 #m, distance from nozzle to skimmer, s in text
s_coll = 59E-3 #m, distance to second skimmer/collimator, d in text

d_mol = 330E-12 #m, diameter of CO2 https://en.wikipedia.org/wiki/Kinetic_diameter
sigma = d_mol #m, collision diameter of molecules, is fitted in paper
gamma = 1.2885 #specific heat ratio cp/cv of CO2

# v_avg = 560 #m/s, mean velocity in the beam, v bar in text

coll_eff = 1 #collision effectiveness

def main():
    v_avg = calc_v_avg(1)
    s_skim = np.arange(4, 30, 0.2)*d_nozzle
    mach = calc_mach(s_skim)
    plt.plot(s_skim/d_nozzle, mach)
    plt.xlabel('Nozzle-skimmer distance (nozzle diameters)')
    plt.ylabel('Mach number')
    plt.axis([0, None, 0, None])
    plt.show()
    plt.close()
    

    density = calc_density(mach, n_0, gamma)
    I_beam = density * v_avg
    plt.plot(s_skim/d_nozzle, I_beam)
    plt.xlabel('Nozzle-skimmer distance (nozzle diameters)')
    plt.ylabel('Beam intensity in second skimmer')
    plt.axis([0, None, 0, None])
    plt.show()
    plt.close()

    s_skim = np.arange(1E-3, 11E-3, 2E-4)
    mach = calc_mach(s_skim)
    density = calc_density(mach, n_0, gamma)
    I_beam = density * v_avg
    plt.plot(s_skim/d_nozzle, I_beam)
    plt.xlabel('Nozzle-skimmer distance (nozzle diameters)')
    plt.ylabel('Beam intensity in second skimmer')
    plt.axis([0, None, 0, None])
    plt.show()
    plt.close()
    print(v_avg)

    s_skim = 10E-3 #m
    co2_fraction = np.arange(0,1.00001,0.01)
    co2_fraction = 1
    n_0_mixed = n_0*co2_fraction
    v_avg = calc_v_avg(co2_fraction) #note: fix the 5/2 here
    gamma_avg = calc_gamma_avg(co2_fraction)
    print('v_avg', v_avg, 'gamma_avg', gamma_avg)

    mach = calc_mach(s_skim)
    density = calc_density(mach, n_0_mixed, gamma_avg)

    I_beam = density * v_avg
    # print(I_beam)

    plt.plot(co2_fraction, I_beam)
    plt.xlabel('CO2 fraction')
    plt.ylabel('Beam intensity in second skimmer')
    plt.axis([0, None, 0, None])
    plt.show()
    plt.close()

    print('pure co2: ', I_beam[-1])



    #Vary skimmer distance

def calc_v_avg(co2_fraction, m1=44, m2=4):
    """
    assumes CO2 (1) in He (2)
    https://www.engineeringtoolbox.com/specific-heat-capacity-gases-d_159.html
    """
    fraction_2 = 1 - co2_fraction
    avg_mass = co2_fraction * m1 + fraction_2 * m2

    E = 5/2*(m1/avg_mass)*R*T_nozzle
    return np.sqrt(2*E/(m1/1000)) #/1000 to get m1 in kg

def calc_gamma_avg(co2_fraction):
    cp1 = 0.844
    cv1 = 0.655
    cp2 = 5.19
    cv2 = 3.12
    fraction_2 = 1 - co2_fraction
    gamma_avg = (co2_fraction * cp1 + fraction_2 * cp2)/(co2_fraction * cv1 + fraction_2 * cv2) #eq 3.5
    return gamma_avg



def calc_mach(s_skim):
    """
    MACH NUMBER from anderson & fenn 1965
    s_skim = number or array containing nozzle-skimmer distances
    pressure_0 = pressure in nozzle
    """
    mfp = kB * T_nozzle / (np.sqrt(2) * np.pi * pressure_0 * d_mol**2) #m, mean free path in nozzle
    kn = mfp/d_nozzle
    mach_1 = 1.23*coll_eff**0.4*kn**(-0.4) #for gamma = 5/3, where dmach/ds = continuum mach rate of change
    s_1 = 0.238*coll_eff**0.6*kn**(-0.6) #for gamma = 5/3


    mach = mach_1 + 0.195*coll_eff/kn/s_1 - 0.195*coll_eff/kn*d_nozzle/s_skim  #mach number at skimmer
    # print('mach_skimmer', mach)

    iscontinuum = s_skim/d_nozzle < s_1
    if isinstance(s_skim, np.ndarray):
        mach[iscontinuum] = 3.2*(s_skim[iscontinuum]/d_nozzle)**(2/3)
    elif iscontinuum: #if it is not an array but in the continuum regime
        mach = 3.2*(s_skim[iscontinuum]/d_nozzle)**(2/3)

    # print('mach fixed', mach)

    mach_t1 = 2.05*coll_eff**0.4*kn**(-0.4)
    # print('mach_t1', mach_t1)
    mach_t = mach_1 + 0.195*coll_eff/kn/s_1
    # print ('mach_t', mach_t, 'mach_1', mach_1)
    print('warning: testing with fractions of mach number')
    return mach

def calc_density(mach, n_0, gamma):
    """
    #for equation 2.14
    mach = number or array of mach numbers
    s_coll = distance to second skimmer
    """

    top = n_0 * mach**2 * d_skim**2 * gamma / (8 * s_coll**2)
    bottom_1 = (1 + (gamma-1)*mach**2/2)**(1/(gamma-1))
    bottom_2 = np.pi * sigma**2 * d_skim * n_0 / np.sqrt(2)
    bottom_3 = 1-np.sqrt(gamma/2/np.pi)*mach*d_skim/2/s_coll
    bottom = bottom_1 + bottom_2*bottom_3
    n_coll = top/bottom #number density at collimator / second skimmer
    # print('number density in second skimmer', n_coll)

    return n_coll



main()

# #-------------------N2-----------------------------------

# T_nozzle = 300 #Kelvin, temperature in nozzle
# d_nozzle = 0.1E-3 #m, nozzle diameter

# pressure_0 = 13.3E2 #pa, pressure in nozzle
# n_0 = pressure_0 / kB / T_nozzle #/m3, number density in nozzle

# d_skim = 1E-3 #m, skimmer diameter, a in text
# s_skim = 3E-3 #m, distance from nozzle to skimmer, s in text
# s_coll = 30E-5 #m, distance to second skimmer/collimator, d in text

# d_mol = 880E-12 #m, diameter of N2 according to paper
# sigma = d_mol #m, collision diameter of molecules, is fitted in paper
# gamma = 1.4 #specific heat ratio cp/cv of N2

# v_avg = 700 #m/s, mean velocity in the beam, v bar in text

# coll_eff = 1 #collision effectiveness

# #-------------------Ar-----------------------------------

# T_nozzle = 300 #Kelvin, temperature in nozzle
# d_nozzle = 0.1E-3 #m, nozzle diameter

# pressure_0 = 1000E2 #pa, pressure in nozzle
# n_0 = pressure_0 / kB / T_nozzle #/m3, number density in nozzle

# d_skim = 1E-3 #m, skimmer diameter, a in text
# s_skim = 3E-3 #m, distance from nozzle to skimmer, s in text
# s_coll = 30E-5 #m, distance to second skimmer/collimator, d in text

# d_mol = 560E-12 #m, diameter of N2 according to paper
# sigma = d_mol #m, collision diameter of molecules, is fitted in paper
# gamma = 1.67 #specific heat ratio cp/cv of N2

# v_avg = 600 #m/s, mean velocity in the beam, v bar in text

# coll_eff = 1 #collision effectiveness


#---------other things

# alpha = np.sqrt(m_gas*v_avg**2/(8*kB*T_beam))*d_skim #some parameter that makes calculation easier
#lambda*n**2/2 is the collision frequency per unit volume at density n

# lamb = 2*sigma**2*np.sqrt(np.pi*R*T_beam/M_gas) #collision parameter

# n_skim = n_0 / (1+(gamma-1)/2*mach**2)**(1/(gamma-1)) #eq 2.10
# n_coll = n_skim*alpha**2/s_coll**2/(1+lamb*n_skim/v_avg*(np.sqrt(np.pi)*alpha-alpha**2/s_coll)) #eq 2.12

# v_trans = -1 #transverse speed of the molecules, w in text

# x = -1 #variable for distance from skimmer

# n_s = -1 #number density at the skimmer

# T_beam = -1 #temperature in the beam, T in text

# c_p = -1 #molar heat capacity


# m_gas = 1 #molecular mass of gas, m in text
# M_gas = 1 #molecular mass of gas, M in text