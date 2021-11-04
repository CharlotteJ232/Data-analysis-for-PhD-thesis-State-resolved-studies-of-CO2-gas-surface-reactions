# To simulate power fluence curves with gaussian velocities etc

import numpy as np
from matplotlib import pyplot as plt
from scipy.integrate import solve_ivp
from tqdm import tqdm 
from multiprocessing import cpu_count
from joblib import Parallel, delayed
import os

# the Optical Bloch equation
def bloch_diffequ_old_unchanged(t,y,E0,w0x,w0y,wzx,wzy,vx,vy,dylaser,dzlaser,mu21,h,z0,vz,zrx,lam,cg,Ecut,lens): #still has cg*mu21, but cg is already included in mu21 in my code
    """
    x=vx*t 
    y=vy*t+dylaser
    z=dzlaser+vz*t+z0

    Rz=(dzlaser+vz*t+z0+(zrx*zrx)/(dzlaser+vz*t+z0))=z+zrx*zrx/z
    wzx=((w0x*np.sqrt(1+(dzlaser+vz*t+z0)*(dzlaser+vz*t+z0)/(zrx*zrx))=w0x*np.sqrt(1+z*z/(zrx/zrx))

    Frequencies in rad/s, and Bloch equations are:
    dy[0]/dt=-(doppler-sweep)*y[1]
    dy[1]/dt=(doppler-sweep)*y[0]+omega(z)*y[2]
    dy[2]/dt=-omega(z)*y[1]fi

    doppler=vz/lam
    sweep=vx*vx*t/(Rz*lam)
    omega(z)=(mu*E/h)*np.sqrt(wx_0*wy_0/(wx_z*wy_z)
    E=E_0*exp(-x*x/(wx_z*wx_z)-y*y/(wy_z*wy_z))

    Laser beam defined as travelling along z axis, molecular beam defined
    as travelling along x axis.
    """
    dydx = np.zeros(3)
    if (E0*(np.exp(-vx*t*vx*t/(wzx*wzx)-((dylaser+vy*t)*(dylaser+vy*t)/(wzy*wzy)))) > Ecut):
        dydx[0]=-2*np.pi*(vz/lam-vx*vx*t/(lam*(dzlaser+vz*t+z0+(zrx*zrx)/(dzlaser+vz*t+z0))))*y[1]
        
        dydx[1]= (2*np.pi*((vz/lam-vx*vx*t/(lam*(dzlaser+vz*t+z0+(zrx*zrx)/(dzlaser+vz*t+z0))))*y[0]+
                (cg*mu21*E0/h)*np.exp(-vx*t*vx*t/((w0x*np.sqrt(1+(dzlaser+vz*t+z0)*(dzlaser+vz*t+z0)/(zrx*zrx)))*
                (w0x*np.sqrt(1+(dzlaser+vz*t+z0)*(dzlaser+vz*t+z0)/(zrx*zrx))))-(dylaser+vy*t)*(dylaser+vy*t)/
                (wzy*wzy))*np.sqrt((w0x*w0y)/((w0x*np.sqrt(1+(dzlaser+vz*t+z0)*(dzlaser+vz*t+z0)/(zrx*zrx)))*wzy))*y[2]))
        
        dydx[2]=(-2*np.pi*((cg*mu21*E0/h)*np.exp(-vx*t*vx*t/((w0x*np.sqrt(1+(dzlaser+vz*t+z0)*(dzlaser+vz*t+z0)/(zrx*zrx)))*
                (w0x*np.sqrt(1+(dzlaser+vz*t+z0)*(dzlaser+vz*t+z0)/(zrx*zrx))))-((dylaser+vy*t)*(dylaser+vy*t)/
                (wzy*wzy)))*np.sqrt((w0x*w0y)/((w0x*np.sqrt(1+(dzlaser+vz*t+z0)*(dzlaser+vz*t+z0)/(zrx*zrx)))*wzy)))*y[1])
    
    return dydx

def bloch_diffequ_old(t,y,E0,w0x,w0y,wzx,wzy,vx,vy,dylaser,dzlaser,mu21,h,z0,vz,zrx,lam,cg,Ecut,lens): #only removed the cg factor in cg*mu21
    """
    x=vx*t 
    y=vy*t+dylaser
    z=dzlaser+vz*t+z0

    Rz=(dzlaser+vz*t+z0+(zrx*zrx)/(dzlaser+vz*t+z0))=z+zrx*zrx/z
    wzx=((w0x*np.sqrt(1+(dzlaser+vz*t+z0)*(dzlaser+vz*t+z0)/(zrx*zrx))=w0x*np.sqrt(1+z*z/(zrx/zrx))

    Frequencies in rad/s, and Bloch equations are:
    dy[0]/dt=-(doppler-sweep)*y[1]
    dy[1]/dt=(doppler-sweep)*y[0]+omega(z)*y[2]
    dy[2]/dt=-omega(z)*y[1]fi

    doppler=vz/lam
    sweep=vx*vx*t/(Rz*lam)
    omega(z)=(mu*E/h)*np.sqrt(wx_0*wy_0/(wx_z*wy_z)
    E=E_0*exp(-x*x/(wx_z*wx_z)-y*y/(wy_z*wy_z))

    Laser beam defined as travelling along z axis, molecular beam defined
    as travelling along x axis.
    """
    dydx = np.zeros(3)
    if (E0*(np.exp(-vx*t*vx*t/(wzx*wzx)-((dylaser+vy*t)*(dylaser+vy*t)/(wzy*wzy)))) > Ecut):
        dydx[0]=-2*np.pi*(vz/lam-vx*vx*t/(lam*(dzlaser+vz*t+z0+(zrx*zrx)/(dzlaser+vz*t+z0))))*y[1]
        
        dydx[1]= (2*np.pi*((vz/lam-vx*vx*t/(lam*(dzlaser+vz*t+z0+(zrx*zrx)/(dzlaser+vz*t+z0))))*y[0]+
                (mu21*E0/h)*np.exp(-vx*t*vx*t/((w0x*np.sqrt(1+(dzlaser+vz*t+z0)*(dzlaser+vz*t+z0)/(zrx*zrx)))*
                (w0x*np.sqrt(1+(dzlaser+vz*t+z0)*(dzlaser+vz*t+z0)/(zrx*zrx))))-(dylaser+vy*t)*(dylaser+vy*t)/
                (wzy*wzy))*np.sqrt((w0x*w0y)/((w0x*np.sqrt(1+(dzlaser+vz*t+z0)*(dzlaser+vz*t+z0)/(zrx*zrx)))*wzy))*y[2]))
        
        dydx[2]=(-2*np.pi*((mu21*E0/h)*np.exp(-vx*t*vx*t/((w0x*np.sqrt(1+(dzlaser+vz*t+z0)*(dzlaser+vz*t+z0)/(zrx*zrx)))*
                (w0x*np.sqrt(1+(dzlaser+vz*t+z0)*(dzlaser+vz*t+z0)/(zrx*zrx))))-((dylaser+vy*t)*(dylaser+vy*t)/
                (wzy*wzy)))*np.sqrt((w0x*w0y)/((w0x*np.sqrt(1+(dzlaser+vz*t+z0)*(dzlaser+vz*t+z0)/(zrx*zrx)))*wzy)))*y[1])
    
    return dydx

def bloch_diffequ(t,y,E0,w0x,w0y,wzx,wzy,vx,vy,dylaser,dzlaser,mu21,h,z0,vz,zrx,lam,cg,Ecut,lens): #refactored from original function
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

def bloch_diffequ_new(t,y,E0,w0x,w0y,wzx,wzy,vx,vy,dylaser,dzlaser,mu21,h,z0,vz,zrx,lam,cg,Ecut,lens): #made from scratch, based on the comments in original function
    x=vx*t 
    y_coordinate=dylaser+vy*t
    z=dzlaser+vz*t+z0

    Rz=z+zrx**2/z
    wzx=w0x*np.sqrt(1+z**2/zrx**2)
    # if not lens:
    #     wzx = w0x

    E=E0*np.exp(-x**2/wzx**2-y_coordinate**2/wzy**2)

    #Only these three are in the final bloch equations
    doppler=vz/lam
    sweep=vx**2*t/(Rz*lam)
    # if not lens:
    #     sweep=0
    omegaz=(mu21*E/h)*np.sqrt(w0x**2/wzx**2)

    dydx = np.zeros(3)
    if E > Ecut:        
        dydx[0]=-2*np.pi*(doppler-sweep)*y[1]
        dydx[1]=2*np.pi*((doppler-sweep)*y[0]+omegaz*y[2])
        dydx[2]=-2*np.pi*omegaz*y[1]

    return dydx



def simulate_trajectory():
    """
    If molecule does not go through skimmers, returns probability -1
    If molecule does go through skimmers, returns calculated probability 
    """
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
            XY = solve_ivp(bloch_diffequ,[-3*T,3*T],[0,0,-1], 
                    args=(E0,w0x,w0y,wzx,wzy,vx,vy,dylaser,dzlaser,mu21,h,z0,vz,zrx,lam,cg,Ecut,lens))
            
            prob=(1+XY['y'][2,:])/2              # Prob=(1+y[2])/2
            return prob[-1]

folderstart = 'C:/Users/jansenc3/Surfdrive/DATA/'
folderstart = 'C:/Users/Werk/surfdrive/DATA/'
savefolder = folderstart+'Laser/Simulations/CO2/'
filenameaddition = '_bloch_new_nolens'
save=False
if save:
    if not os.path.exists(savefolder):
        os.makedirs(savefolder)

## Expt input variables

fx_list=[0.2]       # focus length of lens in xz plane (m)
z0_list=[0.28]  # position of molecular beam from focussed waist (m)
lens = True # false for calculations without lens

vel=2128       # velocity of molecular beam (m/s) (587 for pure)
fwhmv=0.17     # full width half maximum of velocity distribution (as fraction of v)

lam=4255.469316E-9 # wavelength (m) (4252.7)
A21=2/np.pi * 1.407E2  # Einstein coefficient (Hz) (1.811E+02)
cg_list = [1] # transition probability (Clebsch-Gordan coefficient) 0.507,-0.478091,0.377964

ang=0.15       # angular spread of beam out of nozzle (degrees)
rx_list=[0.002]       # radius of unfocused laser in x (m)
ry_list=rx_list      # radius of unfocused laser in y direction (m). coupled to rx_list for the for loop, so there are n iterations and not n^2

#Numbers that change speed and accuracy of the calculation
pcut=5e-05 # minimum laser power for integration (W)
pstep=0.002 # laser power step for fluence curve (W)  
pmin=pstep*0.1     # minimum laser power for fluence curve (W)
pmax=0.05     # maximum laser power for fluence curve (W)

nmax=500  # number of trajectories, was 500

noztolaser=350e-03 # distance between nozzle and laser (m)

#Skimmer data, not up to date, but should probably not be necessary
noztoskim=60e-03   # distance between nozzle and 1st skimmer (m)
noztoskim2=160e-03 # distance between nozzle and 2nd skimmer (m)
znozzle=0 # offset of nozzle to skimmers along laser axis (m)
ynozzle=0 # offset of nozzle to skimmers perpendicular to laser axis (m)
skimdiam=1e-03          # 1st skimmer diameter (m)
skimdiam=skimdiam/2     # 1st skimmer radius (m)
skimdiam2=3e-03         # 2nd skimmer diameter (m)
skimdiam2=skimdiam2/2   # 2nd skimmer radius (m) 

## constants
eps0=8.85e-12	# vacuum permitivity (F/m)
h=6.626e-34    # Plancks constant (Js)
c= 2.997E8        # speed of light (m/s)

#calculate number of cores for parallel processing
num_cores = cpu_count()

## Looping over trajectories
#  for P=pmin:pstep:pmax      # looping over laser power
for rx, ry in zip(rx_list,ry_list):
    for fx in fx_list:
        for z0 in z0_list: 
            for cg in cg_list:
                print('fx ',fx,', z0',z0,', cg ',cg) #print the parameters of this iteration of the loop
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
                mu21=cg*np.sqrt((3*eps0*h*c*c*c*A21)/(2*(wb*wb*wb))) # dipole moment (Cm)
                for P in tqdm(np.arange(pmin, pmax, pstep)):

                    I=2*P/A                  # maximum laser intensity at waist (W/m^2)
                    E0=np.sqrt((2*I)/(eps0*c))  # electic field at waist (V/m)

                    results = Parallel(n_jobs=num_cores)(delayed(simulate_trajectory)() for n in range(nmax))
                    results = np.array(results)
                    # print(results)
                    ispositive = results > 0 
                    probtot = np.average(results[ispositive]) # averaging probability over all trajectories
                        
                    pop.append(probtot)                   # making probability an array to plot later
                    power.append(1000*P)                    # making power an array to plot later
                # end                                    # end of power loop

                plt.plot(power,pop)                 # plots fluence curve
                plt.show()
                plt.close()

                filename='fx_'+str(fx)+'_z0_'+str(z0)+'_lens_'+str(lens)+'_vel_'+str(vel)+'_fwhmv_'+str(fwhmv)+'_lam_'+str(lam)+'_cg_'+str(cg)+'_ang_'+str(ang)+'_rx_'+str(rx)+'_ry_'+str(ry)+'_nmax_'+str(nmax)+'.txt'
                if save:
                    pop = np.array(pop)
                    power = np.array(power)
                    savedata = np.column_stack((power,pop))
                    np.savetxt(savefolder+filename,savedata)
                    # f=open(savefolder+filename,'a')
                    # np.savetxt(f,savedata)
                    # f.close()
print('Done!')


