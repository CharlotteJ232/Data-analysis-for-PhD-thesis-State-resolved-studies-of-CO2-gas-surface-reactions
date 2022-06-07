import numpy as np
from matplotlib import pyplot as plt
import os

folderstart = 'C:/Users/jansenc3/surfdrive/DATA/'
folderstart = 'C:/Users/Werk/surfdrive/DATA/'

folder = folderstart+'Laser/Simulations/CO2/220603/'
savefolder = folder+'/figures/100/'
save=True

transitions = ['R0', 'R2', 'R4', 'R6', 'R8', 'R10']
# transitions = ['R4']

cg_lists = {'R0': [1],
            'R2': [0.507,0.478091,0.377964],
            'R4': [0.389249, 0.381385, 0.356753, 0.3114, 0.23355],
            'R6': [0.328165, 0.324799, 0.314485, 0.229668, 0.269309, 0.229668, 0.169031],
            'R8': [0.289122, 0.287331, 0.281892, 0.272587, 0.258997, 0.240399, 0.215499, 0.181724, 0.132453],
            'R10':[0.261387, 0.260304, 0.25703, 0.251478, 0.243492, 0.232823, 0.219079, 0.201631, 0.179402, 0.150287, 0.108893]}

# cg_lists = {'Test':[0.75],
#             'Test_avg':[0.5, 1],
#             'R0': [1],
#             'R2': [0.774597 ,0.730297,0.57735],
#             'R4': ['']}

wl = {  'Test': 4.255469316e-06,
        'Test_avg': 4.255469316e-06,
        'R0': 4.255469316e-06,
        'R2': 4.2527e-06,
        'R4': 4.249979e-06,
        'R6': 4.2473e-06,
        'R8': 4.24468e-06,
        'R10':4.2421e-06}

vel=587       # velocity of molecular beam (m/s)
fwhmv=0.17      # full width half maximum of velocity distribution (as fraction of v)
ang=0.15       # angular spread of beam out of nozzle (degrees)
# ang=0.3437705518714731
rx=3e-03       # radius of unfocused laser in x (m)
ry=3e-03       # radius of unfocused laser in y direction (m)
fx=0.2
z0=0.28
lens=True
nmax=100  # number of trajectories, was 500

def main():
    if save:
        if not os.path.exists(savefolder):
            os.makedirs(savefolder)

    populations = {}
    for transition in transitions: 
        loadfolder = folder+transition+'/'
        loadfolder = folder
        populations['power'], populations[transition] = average_cg(transition,cg_lists[transition],title=transition, save=save, loadfolder=loadfolder, savefolder=savefolder, lam=wl[transition])



def average_cg_oldweights(cg_list, title='', save=False, loadfolder='', savefolder='', lam=None):
    totalpopulation=[]
    weight = 1 #for the first entry only, corresponding to m1=m2=0
    for cg in cg_list:
        filename='fx_'+str(fx)+'_z0_'+str(z0)+'_lens_'+str(lens)+'_vel_'+str(vel)+'_fwhmv_'+str(fwhmv)+'_lam_'+str(lam)+'_cg_'+str(cg)+'_ang_'+str(ang)+'_rx_'+str(rx)+'_ry_'+str(ry)+'_nmax_'+str(nmax)+'.txt'        
        power, population = np.loadtxt(loadfolder+filename,unpack=True)
        # plt.plot(power,population,label=str(cg))
        totalpopulation.append(weight*np.array(population))
        weight = 2
    totalpopulation = np.sum(np.array(totalpopulation),0)
    totalpopulation /= (2*len(cg_list)-1)
    plt.plot(power, totalpopulation,label='averaged')
    plt.title(title)
    plt.xlabel('Laser power (mW)')
    plt.ylabel('Excited population')
    plt.axis([0,None,0,None])
    if save:
        plt.savefig(savefolder+title+'.png')
        np.savetxt(loadfolder+'averaged.txt', np.column_stack((power, totalpopulation)))
    # plt.legend()
    plt.show()
    plt.close()

    return power, totalpopulation

def average_cg(transition, cg_list, title='', save=False, loadfolder='', savefolder='', lam=None):
    totalpopulation=[]
    for cg in cg_list:
        filename='fx_'+str(fx)+'_z0_'+str(z0)+'_lens_'+str(lens)+'_vel_'+str(vel)+'_fwhmv_'+str(fwhmv)+'_lam_'+str(lam)+'_cg_'+str(cg)+'_ang_'+str(ang)+'_rx_'+str(rx)+'_ry_'+str(ry)+'_nmax_'+str(nmax)+'.txt'
        filename=transition+'_fx_'+str(fx)+'_z0_'+str(z0)+'_lens_'+str(lens)+'_vel_'+str(vel)+'_fwhmv_'+str(fwhmv)+'_lam_'+str(lam)+'_ang_'+str(ang)+'_rx_'+str(rx)+'_ry_'+str(ry)+'_nmax_'+str(nmax)+'.txt'  
        power, population = np.loadtxt(loadfolder+filename,unpack=True)
        # plt.plot(power,population,label=str(cg))
        totalpopulation.append(np.array(population))
    totalpopulation = np.sum(np.array(totalpopulation),0)
    totalpopulation /= len(cg_list)
    plt.plot(power, totalpopulation,label='averaged')
    plt.title(title)
    plt.xlabel('Laser power (mW)')
    plt.ylabel('Excited population')
    plt.axis([0,None,0,None])
    if save:
        plt.savefig(savefolder+title+'.png')
        np.savetxt(loadfolder+transition+'_averaged.txt', np.column_stack((power, totalpopulation)))
    # plt.legend()
    plt.show()
    plt.close()

    return power, totalpopulation


main()