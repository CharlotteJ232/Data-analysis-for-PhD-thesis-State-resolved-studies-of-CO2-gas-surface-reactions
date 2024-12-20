import numpy as np
from matplotlib import pyplot as plt
import os

folderstart = 'C:/Users/jansenc3/surfdrive/DATA/'
folderstart = 'C:/Users/Werk/surfdrive/DATA/'
loadfolder = folderstart+'Laser/Simulations/CO2/2022-03/'

savefolder = folderstart+'Laser/Simulations/CO2/figures/'

transitions = ['R0', 'R2', 'R4']

vel=587       # velocity of molecular beam (m/s)
fwhmv=0.17      # full width half maximum of velocity distribution (as fraction of v)

lam=4252.7E-9 # wavelength (m)
A21=1.811E+02  # Einstein coefficient (Hz)
cg_list = [0.507,-0.478091,0.377964]

ang=0.15       # angular spread of beam out of nozzle (degrees)
# ang=0.3437705518714731
rx=3e-03       # radius of unfocused laser in x (m)
ry=3e-03       # radius of unfocused laser in y direction (m)

fx=0.2
z0=0.28
lens=True

nmax=500  # number of trajectories, was 500

def main():
    # fx_list=[0.1,0.2,0.3,0.4,0.5,0.8]       # focus length of lens in xz plane (m)
    # z0_list=[0.05,0.1,0.28,0.5]  # position of molecular beam from focussed waist (m)
    # lens=True
    # plot_fx_z0(fx_list,z0_list,lens=lens,title='With lens')

    # lens=False
    # fx_list=[5,10,20,50,100,200,1000]
    # z0_list=[0]+fx_list
    # plot_fx_z0(fx_list,z0_list,lens=lens, title='Without lens')

    # lens=True
    # fx=0.2
    # z0=0.28
    # rx_list=[0.003] 
    # ry_list=rx_list
    # plot_r(rx_list,ry_list, fx, z0, lens=lens, title='')

    average_cg(cg_list,title='PED', save=False)


def plot_fx_z0(fx_list,z0_list, average_cg_coefficient=True, lens=True, title=''):

    for fx,c in zip(fx_list,['r', 'b', 'black', 'gray', 'g', 'cyan', 'gold', 'orange','magenta','lightgreen']):
        for z0,ls in zip(z0_list,['--',':','-','-.','-','--',':','-','-.']):
            first = True #for adding different datasets with different cg
            incomplete = False #becomes true when a cg value is missing
            for cg,cg_weight in zip(cg_list,cg_weight_list):
                filename='fx_'+str(fx)+'_z0_'+str(z0)+'_lens_'+str(lens)+'_vel_'+str(vel)+'_fwhmv_'+str(fwhmv)+'_lam_'+str(lam)+'_cg_'+str(cg)+'_ang_'+str(ang)+'_rx_'+str(rx)+'_ry_'+str(ry)+'_nmax_'+str(nmax)+'.txt'
                if os.path.exists(loadfolder+filename):
                    power,population=np.loadtxt(loadfolder+filename,unpack=True)
                    if first:
                        totalpopulation = cg_weight*np.array(population)
                    else: 
                        totalpopulation += cg_weight*np.array(population)
                    first = False
                    if not average_cg_coefficient:
                        plt.plot(power,population,ls=ls, c=c,label='fx='+str(fx)+', z0='+str(z0)+', cg='+str(cg))
                else:
                    incomplete=True
            if average_cg_coefficient and not incomplete:
                totalpopulation /= np.sum(cg_weight_list)
                plt.plot(power,totalpopulation,ls=ls, c=c,label='fx='+str(fx)+', z0='+str(z0))
    plt.legend()
    plt.xlabel('Laser power (mW)')
    plt.ylabel('Excited population')
    plt.title(title)
    plt.axis([0,None,0,None])
    plt.show()
    plt.close()



def plot_r(rx_list,ry_list, fx,z0, average_cg_coefficient=True, lens=True, title=''):
    for rx,ry in zip(rx_list,ry_list):
        first = True #for adding different datasets with different cg
        incomplete = False #becomes true when a cg value is missing
        for cg,cg_weight in zip(cg_list,cg_weight_list):
            filename='fx_'+str(fx)+'_z0_'+str(z0)+'_lens_'+str(lens)+'_vel_'+str(vel)+'_fwhmv_'+str(fwhmv)+'_lam_'+str(lam)+'_cg_'+str(cg)+'_ang_'+str(ang)+'_rx_'+str(rx)+'_ry_'+str(ry)+'_nmax_'+str(nmax)+'.txt'
            if os.path.exists(loadfolder+filename):
                power,population=np.loadtxt(loadfolder+filename,unpack=True)
                if first:
                    totalpopulation = cg_weight*np.array(population)
                else: 
                    totalpopulation += cg_weight*np.array(population)
                first = False
                if not average_cg_coefficient:
                    plt.plot(power,population, label=', rx='+str(rx)+', ry='+str(ry))
            else:
                incomplete=True
        if average_cg_coefficient and not incomplete:
            totalpopulation /= np.sum(cg_weight_list)
            plt.plot(power,totalpopulation,label=', rx='+str(rx)+', ry='+str(ry))
    plt.legend()
    plt.xlabel('Laser power (mW)')
    plt.ylabel('Excited population')
    plt.title(title)
    plt.axis([0,None,0,None])
    plt.show()
    plt.close()



def average_cg(cg_list, title='', save=False):
    ## Expt input variables
    # fx=0.2       # focus length of lens in xz plane (m)
    # z0=0   # position of molecular beam from focussed waist (m)
    # lens = False # false for calculations without lens

    # vel=587       # velocity of molecular beam (m/s)
    # fwhmv=0.25      # full width half maximum of velocity distribution (as fraction of v)

    # lam=4252.7E-9 # wavelength (m)
    # cg = 0.507 # transition probability (Clebsch-Gordan coefficient)
    # cg = -0.478091 #for m1=1
    # # cg = 0.377964 #for m1=2

    # ang=0.15       # angular spread of beam out of nozzle (degrees)
    # rx=2e-03       # radius of unfocused laser in x (m)
    # ry=2e-03       # radius of unfocused laser in y direction (m)

    # nmax=500  # number of trajectories, was 500

    totalpopulation=[]
    weight = 1 #for the first entry only, corresponding to m1=m2=0
    for cg in cg_list:
        filename='fx_'+str(fx)+'_z0_'+str(z0)+'_lens_'+str(lens)+'_vel_'+str(vel)+'_fwhmv_'+str(fwhmv)+'_lam_'+str(lam)+'_cg_'+str(cg)+'_ang_'+str(ang)+'_rx_'+str(rx)+'_ry_'+str(ry)+'_nmax_'+str(nmax)+'.txt'
        power, population = np.loadtxt(loadfolder+filename,unpack=True)
        # plt.plot(power,population,label=str(cg))
        totalpopulation.append(weight*np.array(population))
        weight = 2
    totalpopulation = np.sum(np.array(totalpopulation),0)
    plt.plot(power, totalpopulation/(2*len(cg_list)-1),label='averaged')
    plt.title(title)
    plt.xlabel('Laser power (mW)')
    plt.ylabel('Excited population')
    plt.axis([0,None,0,None])
    if save:
        plt.savefig(loadfolder+title+'.png')
    # plt.legend()
    plt.show()
    plt.close()

    return 



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