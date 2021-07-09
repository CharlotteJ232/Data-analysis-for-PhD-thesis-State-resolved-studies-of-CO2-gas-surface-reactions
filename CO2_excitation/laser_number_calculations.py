import numpy as np
    
z = np.array([0, 28, 43, 58])*1E-2 #positions relative to the focal point, in m


def main():
    calculate_laser_diameter(z)

def calculate_laser_diameter(z):
    print('calculating laser diameter..')
    wl = 4252.7E-9 #wavelength in m
    r_laser = 2E-3 #radius of unfocused laser, in m
    f_lens = 20E-2 #focal length of lens in m

    w_focalpoint = f_lens*wl/(r_laser*np.pi)
    rayleigh_range = np.pi*w_focalpoint**2/(wl)
    w_z = w_focalpoint*np.sqrt(1+z**2/rayleigh_range**2)

    for i in range(len(z)):
        print('z=',np.round(z[i]*100,3),'cm, w_z=',np.round(w_z[i]*1000,3),'mm')


def calculate_polarisation():
    n_gold = 2.7022 #at 4252.7 nm wavelength
    n_aluminium = 7.3594

main()

