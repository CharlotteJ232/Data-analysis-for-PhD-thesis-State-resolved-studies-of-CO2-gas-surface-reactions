# To simulate power fluence curves with gaussian velocities etc

import numpy as np
from matplotlib import pyplot as plt
from scipy.integrate import solve_ivp
from tqdm import tqdm 
from multiprocessing import cpu_count
from joblib import Parallel, delayed
import os
import random

# the Optical Bloch equation
def bloch_diffequ(t,y,E0,w0x,w0y,wzx,wzy,vx,vy,dylaser,dzlaser,mu21,h,z0,vz,zrx,lam,Ecut,lens): #refactored from original function
    """
    Laser beam defined as travelling along z axis, molecular beam defined
    as travelling along x axis.
    """
    x=vx*t
    y_coordinate=dylaser+vy*t
    z=dzlaser+vz*t+z0

    Rz=z+zrx*zrx/z
    if lens:
        wzx=w0x*np.sqrt(1+z**2/zrx**2)
    else:
        wzx=w0x

    E=E0*np.exp(-x**2/wzx**2-y_coordinate**2/wzy**2)

    doppler=vz/lam
    if lens:
        sweep=vx**2*t/(Rz*lam)
    else:
        sweep=0
    omegaz=(mu21*E/h)*np.sqrt(w0x*w0y/(wzx*wzy))

    dydx = np.zeros(3)
    if (E0*(np.exp(-x*x/(wzx*wzx)-((y_coordinate)*(y_coordinate)/(wzy*wzy)))) > Ecut):
        dydx[0]=-2*np.pi*(doppler-sweep)*y[1] 
        dydx[1]= 2*np.pi*((doppler-sweep)*y[0]+omegaz*y[2])
        dydx[2]=-2*np.pi*omegaz*y[1]
    
    return dydx

def random_CG(j):
    """
    returns cg for random m
    < j1 m 1 0 | j1+1 m >
    https://en.wikipedia.org/wiki/Clebsch%E2%80%93Gordan_coefficients#Special_cases
    """
    m_array = np.arange(-j, j+0.1, 1)
    m = random.choice(m_array)
    upper = (j-m+1)*(j+m+1)
    lower = (2*j+1)*(j+1)
    return(np.sqrt(upper/lower))
    
def random_mu_direction():
    """
    Returns the projection of a unit vector with a random orientation to the z axis. 
    Apparently the distribution of projections on the z axis is uniform
    """
    return np.random.random()

def random_mu_rotation():
    """
    projection of randomly oriented circle on the polarization axis
    """
    theta = np.arccos(np.random.random())
    avg = 2/np.pi * np.sin(theta)
    return avg


def simulate_trajectory():
    """
    If molecule does not go through skimmers, returns probability -1
    If molecule does go through skimmers, returns calculated probability 
    """
    #cg random selection
    cg = random_CG(j)
    if j == 0:
        factor = random_mu_direction()
    else:
        factor = random_mu_rotation()
    factor = 1
    mu21=cg*factor*np.sqrt((3*eps0*h*c*c*c*A21)/(2*(wb*wb*wb))) # dipole moment (Cm)

    # velocity selection
    sintheta=np.random.rand(1)       # polar angle of v - need to sample from sin distribution due to sampling over sphere
    theta=np.arcsin(sintheta)       # in radians
    theta=theta*ang/90              # 
    phi=np.random.rand(1)*2*np.pi         # azimuthal angle of v, in radians
    v=np.random.randn(1)*fwhmv*vel+vel # velocity
    vz=v*np.sin(theta)*np.cos(phi) # initial velocity in z direction
    vy=v*np.sin(theta)*np.sin(phi) # initial velocity in y direction
    vx=v*np.cos(theta)          # initial velocity in x direction

    # 1st skimmer
    tskim=noztoskim/vx                         # time to reach 1st skimmer
    dyskim=tskim*vy-ynozzle                    # y position wrt 1st skimmer
    dzskim=tskim*vz-znozzle                    # z position wrt 1st skimmer
    drskim=np.sqrt(dyskim*dyskim+dzskim*dzskim)   # determining if trajectory goes through 1st skimmer

    # 2nd skimmer
    tskim2=noztoskim2/vx                       # time to reach 2nd skimmer
    dyskim2=tskim2*vy-ynozzle                  # y position wrt 2nd skimmer
    dzskim2=tskim2*vz-znozzle                  # z position wrt 2nd skimmer
    drskim2=np.sqrt(dyskim2*dyskim2+dzskim2*dzskim2) # determining if trajectory goes through 2nd skimmer

    if (drskim > skimdiam):
        iskim=iskim+1            # molecule doesn't go through 1st skimmer
        return -1
    else:
        if (drskim2 > skimdiam2):  
            iskim=iskim+1          # molecule doesn't go through 2nd skimmer
            return -1
        else:
            tlaser=noztolaser/vx           # time taken to reach laser centre (s)
            dylaser=tlaser*vy              # distance of molecule perpendicular to plane of laser (m)
            dzlaser=tlaser*vz              # displacement of molecule from z0 (m)
            z=z0+dzlaser                   # distance of molecule from focussed waist (m)
            wzx=rx #when not using lens
            if lens:
                wzx=w0x*np.sqrt(1+z*z/(zrx*zrx))  # beam waist at position z (m)
            wzy=ry                         # beam not focussed in y direction (m)
            T=wzx/vx                       # half transit time of beam
            Icut=2*pcut/(wzx*wzy)          # intensity cut off for integration
            Ecut=np.sqrt((2*Icut)/(eps0*c))   # electric field cut off for integration (approx 200 works well)

            ## Solving the Bloch equations
            XY = solve_ivp(bloch_diffequ,[-3*T,3*T],[0,0,-1], 
                    args=(E0,w0x,w0y,wzx,wzy,vx,vy,dylaser,dzlaser,mu21,h,z0,vz,zrx,lam,Ecut,lens))
            
            prob=(1+XY['y'][2,:])/2              # Prob=(1+y[2])/2
            return prob[-1]



transitions = ['R0', 'R2', 'R4', 'R6', 'R8', 'R10']
# transitions = ['R2']


wl = {  'Test': 4.255469316e-06,
        'R0': 4.255469316e-06,
        'R2': 4.2527e-06,
        'R4': 4.249979e-06,
        'R6': 4.2473e-06,
        'R8': 4.24468e-06,
        'R10':4.2421e-06}

A21_dic = { 'Test': 140.70,
            'R0': 140.70,
            'R2': 181.10,
            'R4': 192.40,
            'R6': 197.70,
            'R8': 201.00,
            'R10':203.10}


folderstart = 'C:/Users/jansenc3/Surfdrive/DATA/'
folderstart = 'C:/Users/Werk/surfdrive/DATA/'
folder = folderstart+'Laser/Simulations/CO2/231206/test/'

filenameaddition = ''
save=True
save_figures = True


## Expt input variables
fx_list=[0.2]       # focus length of lens in xz plane (m)
z0_list=[0.28]  # position of molecular beam from focussed waist (m)
lens = True # false for calculations without lens

vel=587      # velocity of molecular beam (m/s) (587 for pure, 954 for 7/12,  2041 for 1/14)
fwhmv_list=[0.17]     # full width half maximum of velocity distribution (as fraction of v)

#new setup with laser after skimmer
noztolaser=20e-03 # distance between nozzle and laser (m)
#old setup with laser before crystal
noztolaser=350e-03 # distance between nozzle and laser (m)
noztoslidingvalve=300e-3 #distande between nozzle and sliding valve in m
holediam=1.8e-3 #diameter of the hole in the sliding valve in m
ang=np.arctan(holediam/noztoslidingvalve)*180/np.pi #for hole in sliding valve (do not use for PED)
ang=0.15       # angular spread of beam out of nozzle (degrees) (0.15 for PED)
rx_list=[0.003, 0.006]       # radius of unfocused laser in x (m)
ry_list=rx_list      # radius of unfocused laser in y direction (m). coupled to rx_list for the for loop, so there are n iterations and not n^2

#Numbers that change speed and accuracy of the calculation
pcut=5e-05 # minimum laser power for integration (W)
pstep=0.002 # laser power step for fluence curve (W)  
pmin=pstep*0.1     # minimum laser power for fluence curve (W)
pmax=0.04     # maximum laser power for fluence curve (W)

nmax=100  # number of trajectories for each laser power step, 1000 for paper

## constants
eps0=8.85e-12	# vacuum permitivity (F/m)
h=6.626e-34    # Plancks constant (Js)
c= 2.997E8        # speed of light (m/s)

#Skimmer data, not up to date, but should probably not be necessary
noztoskim=60e-03   # distance between nozzle and 1st skimmer (m)
noztoskim2=160e-03 # distance between nozzle and 2nd skimmer (m)
znozzle=0 # offset of nozzle to skimmers along laser axis (m)
ynozzle=0 # offset of nozzle to skimmers perpendicular to laser axis (m)
skimdiam=1e-03          # 1st skimmer diameter (m)
skimdiam=skimdiam/2     # 1st skimmer radius (m)
skimdiam2=3e-03         # 2nd skimmer diameter (m)
skimdiam2=skimdiam2/2   # 2nd skimmer radius (m) 

#calculate number of cores for parallel processing
num_cores = cpu_count()-1 #-1 to save one core for other tasks

##START OF SCRIPT## 

##looping over different transitions
for transition in transitions:
    print(transition)
    j = int(transition[1:])
    lam = wl[transition]
    A21 = A21_dic[transition]
    savefolder = folder
    if save:
        if not os.path.exists(savefolder):
            os.makedirs(savefolder)

    for rx, ry in zip(rx_list,ry_list):
        for fx in fx_list:
            for z0 in z0_list: 
                for fwhmv in fwhmv_list:
                    print('fx ',fx,', z0',z0) #print the parameters of this iteration of the loop
                    pop = [] #make empty population list
                    power = [] #make empty power list
                    ## calculating values that are the same for each trajectory
                    w0x=rx #for no lens
                    if lens:
                        w0x=lam*fx/(np.pi*rx)   # x beam waist at focal point of laser (m)
                    w0y=ry                  # laser beam not focussed in y direction
                    zrx=np.pi*w0x*w0x/lam   # Rayleigh range before focussing (m)
                    A=np.pi*w0x*w0y            # area of focussed laser beam at waist (m^2)
                    wb=2*np.pi*c/lam        # transition frequency

                    ## Looping over laser power
                    for P in tqdm(np.arange(pmin, pmax, pstep)):
                        I=2*P/A                  # maximum laser intensity at waist (W/m^2)
                        E0=np.sqrt((2*I)/(eps0*c))  # electic field at waist (V/m)

                        ## Looping over trajectories
                        results = Parallel(n_jobs=num_cores)(delayed(simulate_trajectory)() for n in range(nmax))
                        results = np.array(results)
                        # print(results)
                        ispositive = results > 0 
                        probtot = np.average(results[ispositive]) # averaging probability over all trajectories
                            
                        pop.append(probtot)                   # making probability an array to plot later
                        power.append(1000*P)                    # making power an array to plot later
                    # end                                    # end of power loop

                    filename=transition+'_fx_'+str(fx)+'_z0_'+str(z0)+'_lens_'+str(lens)+'_vel_'+str(vel)+'_fwhmv_'+str(fwhmv)+'_lam_'+str(lam)+'_ang_'+str(ang)+'_rx_'+str(rx)+'_ry_'+str(ry)+'_nmax_'+str(nmax)


                    plt.plot(power,pop)                 # plots fluence curve
                    plt.axis([0, pmax*1000, 0, 1])
                    plt.xlabel('Laser power (mW)')
                    plt.ylabel('Excited population')
                    if save_figures:
                        plt.savefig(savefolder+filename+'.png')

                    plt.show()
                    plt.close()

                    if save:
                        pop = np.array(pop)
                        power = np.array(power)
                        savedata = np.column_stack((power,pop))
                        np.savetxt(savefolder+filename+'.txt',savedata)
                        # f=open(savefolder+filename,'a')
                        # np.savetxt(f,savedata)
                        # f.close()
print('Done!')

