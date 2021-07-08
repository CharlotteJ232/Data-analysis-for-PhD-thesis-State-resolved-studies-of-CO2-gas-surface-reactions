import numpy as np
from matplotlib import pyplot as plt

datadic = {0.7:{'kw_day':'03 Mar/200310', 'TOF_day':'03 Mar/200310','c':'red'},
           0.8:{'kw_day':'05 May/200512', 'TOF_day':'01 Jan/200128','c':'orange'},
           1.2:{'kw_day':'03 Mar/200305', 'TOF_day':'02 Feb/200220','c':'gold'},
           2.0:{'kw_day':'06 Jun/200602', 'TOF_day':'05 May/200511','c':'yellow'},
           2.9:{'kw_day':'03 Mar/200303', 'TOF_day':'02 Feb/200218','c':'greenyellow'},
           4.0:{'kw_day':'03 Mar/200311', 'TOF_day':'03 Mar/200306','c':'green'},
           4.77:{'kw_day':'05 May/200519', 'TOF_day':'05 May/200519','c':'mediumseagreen'},
           6.1:{'kw_day':'05 May/200520', 'TOF_day':'02 Feb/200204','c':'teal'},
           7.8:{'kw_day':'08 Aug/200817', 'TOF_day':'08 Aug/200813_2','c':'blue'},
           9.4:{'kw_day':'08 Aug/200814', 'TOF_day':'08 Aug/200813','c':'darkblue'},
           10.7:{'kw_day':'07 Jul/200727', 'TOF_day':'07 Jul/200727','c':'indigo'},
           12.9:{'kw_day':'09 Sep/200904', 'TOF_day':'09 Sep/200904','c':'purple'},
           13.1:{'kw_day':'09 Sep/200910', 'TOF_day':'09 Sep/200910','c':'black'}}

h = 6.626E-34 #planck constant
m = 4 * 1.66E-27 #D2 mass in kg
Ewell = 12 /6.022E23 *1000 #J, attractive potential of the surface
E = Ewell + np.array(list(datadic.keys())) / 6.022E23 * 1000 #from kJ/mol to J, energy of the D2 molecules
p = np.sqrt(2*m*E) #momentum of the D2 molecules
l = h/p #wavelength of the D2 molecules in m
d = 1E-9*np.absolute(1/np.loadtxt("C:/Users/Werk/surfdrive/DATA/D2 sticking probability as a function of step density and kinetic energy/Images/step_density.txt",usecols=2,skiprows=1))
#step density ^
lmesh, dmesh = np.meshgrid(l,d)
Emesh, dmesh2 = np.meshgrid(E,d)

print(l)
print(d)

m_list = np.arange(1,15,2) #diffraction peak numbers

for m in m_list:
    theta_m = np.arcsin(m*lmesh/dmesh)

    E_perp = np.cos(theta_m) *Emesh
    trapped = E_perp < Ewell

    print(lmesh[trapped],dmesh[trapped])

    fig = plt.figure()
    ax = plt.axes(projection='3d')
    ax.contour3D(lmesh, dmesh, theta_m, 50, cmap='binary')
    ax.set_xlabel('l')
    ax.set_ylabel('d')
    ax.set_zlabel('theta_m')



