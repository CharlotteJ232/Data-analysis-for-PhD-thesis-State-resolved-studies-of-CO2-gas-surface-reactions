import numpy as np
from matplotlib import pyplot as plt

B = 0.39 #cm-1 or 11.70 GHz = 1.17E10 Hz https:B //cccbdb.nist.gov/exp2x.asp?casno=124389&charge=0 
B = 39 #m-1 #B bar on wikipedia
hc = 1.986E-25 #Jm
c = 3E8
h = hc / c
hbar = 1.05457E-34 #Js
E_0 = B*hc #J #B on wikipedia
I = hbar**2/(2*E_0)
I2 = h / (8 * np.pi**2 * B * c)


def main():
    J = np.array([2,4,6,8,10])
    E = E_0*J*(J+1)
    f = get_period(E, I)
    print(I)
    print(I2)
    print('rotational frequency ',f)

    print("time spent in field", 0.002/587, ' s')


def get_period(E, I):
    """
    E in J
    I in Js^2
    resulting in f in 1/s
    """
    T = np.sqrt(I/(8*np.pi*np.pi*E))
    return T


main()
