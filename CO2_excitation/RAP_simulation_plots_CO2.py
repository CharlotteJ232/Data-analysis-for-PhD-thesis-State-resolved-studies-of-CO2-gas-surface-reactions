import numpy as np
from matplotlib import pyplot as plt

folderstart = 'P:/Surfdrive/DATA/'
folderstart = 'C:/Users/Werk/surfdrive/DATA/'
loadfolder = folderstart+'Laser/Simulations/CO2/'

def main():
    average_cg()

def average_cg():
    ## Expt input variables
    fx=0.2       # focus length of lens in xz plane (m)
    z0=0   # position of molecular beam from focussed waist (m)
    lens = False # false for calculations without lens

    vel=587       # velocity of molecular beam (m/s)
    fwhmv=0.25      # full width half maximum of velocity distribution (as fraction of v)

    lam=4252.7E-9 # wavelength (m)
    cg = 0.507 # transition probability (Clebsch-Gordan coefficient)
    cg = -0.478091 #for m1=1
    # cg = 0.377964 #for m1=2

    ang=0.15       # angular spread of beam out of nozzle (degrees)
    rx=2e-03       # radius of unfocused laser in x (m)
    ry=2e-03       # radius of unfocused laser in y direction (m)

    nmax=500  # number of trajectories, was 500

    totalpopulation=[]
    for cg,weight in zip([0.507,-0.478091,0.377964],[1,2,2]):
        filename='fx_'+str(fx)+'_z0_'+str(z0)+'_lens_'+str(lens)+'_vel_'+str(vel)+'_fwhmv_'+str(fwhmv)+'_lam_'+str(lam)+'_cg_'+str(cg)+'_ang_'+str(ang)+'_rx_'+str(rx)+'_ry_'+str(ry)+'_nmax_'+str(nmax)+'.txt'
        power, population = np.loadtxt(loadfolder+filename,unpack=True)
        plt.plot(power,population,label=str(cg))
        totalpopulation.append(weight*np.array(population))
    totalpopulation = np.array(totalpopulation[0])+np.array(totalpopulation[1])+np.array(totalpopulation[2])
    plt.plot(power, totalpopulation/5,label='averaged')
    plt.axis([0,None,0,None])
    plt.legend()
    plt.show()
    plt.close()



def plot_old():
    f_z0_list = [(0.2, 0.28, 'b', 'lens'), (0.2, 0, 'r', 'no lens')]

    datadic = {}
    for f_z0 in f_z0_list:
        datadic[f_z0] = {}
        datadic[f_z0]['power'], datadic[f_z0]['population'] = np.loadtxt(loadfolder+'f_'+str(f_z0[0])+'_z_'+str(f_z0[1])+'.txt',skiprows=1, unpack=True)
        

    for f_z0_sublist in [f_z0_list]:
        for f_z0 in f_z0_sublist:   
            plt.plot(datadic[f_z0]['power'], datadic[f_z0]['population'], label=str(f_z0[3]), c=f_z0[2])
        plt.axis([0,50,0,1])
        plt.legend()
        plt.xlabel('Laser power (mW)')
        plt.ylabel('Excited population')
        plt.show()
        plt.close()

main()