# To simulate power fluence curves with gaussian velocities etc

import numpy as np
from matplotlib import pyplot as plt
from scipy.integrate import solve_ivp
from tqdm import tqdm 
import os

# the Optical Bloch equation
def bloch_diffequ(t,y,E0,w0x,w0y,wzx,wzy,vx,vy,dylaser,dzlaser,mu21,h,z0,vz,zrx,lam,cg,Ecut):
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



folderstart = 'C:/Users/jansenc3/Surfdrive/DATA/'
# folderstart = 'C:/Users/Werk/surfdrive/DATA/'
savefolder = folderstart+'Laser/Simulations/CO2/'
save=False
if save:
    if not os.path.exists(savefolder):
        os.makedirs(savefolder)

## Expt input variables
lam=4252.7E-9 # wavelength (m)
fx=0.2       # focus length of lens in xz plane (m)
lens = True # false for calculations without lens
rx=2e-03       # radius of unfocused laser in x (m)
ry=2e-03       # radius of unfocused laser in y direction (m)
vel=587       # velocity of molecular beam (m/s)
fwhmv=0.2      # full width half maximum of velocity distribution (as fraction of v)
ang=0.15       # angular spread of beam out of nozzle (degrees)

z0=0   # position of molecular beam from focussed waist (m)
A21=1.811E+02  # Einstein coefficient (Hz)
cg=0.507 # transition probability (Clebsch-Gordan coefficient)

pcut=5e-05 # minimum laser power for integration (W)
pmin=0     # minimum laser power for fluence curve (W)
pmax=0.05     # maximum laser power for fluence curve (W)
pstep=0.001 # laser power step for fluence curve (W)  , was 0.02

nmax=5  # number of trajectories, was 500

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
ip=0           # counter for powers

## calculating values that are the same for each trajectory
w0x=rx #for no lens
if lens:
    w0x=lam*fx/(np.pi*rx)   # x beam waist at focal point of laser (m)
w0y=ry                  # laser beam not focussed in y direction
zrx=np.pi*w0x*w0x/lam   # Rayleigh range before focussing (m)
A=np.pi*w0x*w0y            # area of focussed laser beam at waist (m^2)
wb=2*np.pi*c/lam        # transition frequency
mu21=cg*np.sqrt((3*eps0*h*c*c*c*A21)/(2*(wb*wb*wb))) # dipole moment (Cm)

## Looping over trajectories
pop = []
power = []
#  for P=pmin:pstep:pmax      # looping over laser power
for P in tqdm(np.arange(pmin, pmax, pstep)):
    ip=ip+1                   # Counter for powers used in a later loop
    I=2*P/A                  # maximum laser intensity at waist (W/m^2)
    E0=np.sqrt((2*I)/(eps0*c))  # electic field at waist (V/m)
    iskim=0                  # zeroing counters for molecules that don't make it through the skimmers
    icount=0                 # zeroing counters for molecules that do make it through the skimmer
    probtot=0                # zeroing counters for probability

    # for ntraj=1:1:nmax       # looping over trajectories
    for ntraj in np.arange(nmax):
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
        else:
            if (drskim2 > skimdiam2):  
                iskim=iskim+1          # molecule doesn't go through 2nd skimmer
            else:
                icount=icount+1                # molecule does go through skimmer
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
                # [X,Y] = ode45(@(t,y) bloch_diffequ(t,y,E0,w0x,w0y,wzx,wzy,vx,vy,dylaser,dzlaser,mu21,h,z0,vz,zrx,lam,cg,Ecut),[-3*T,3*T],[0,0,-1])
                XY = solve_ivp(bloch_diffequ,[-3*T,3*T],[0,0,-1], 
                     args=(E0,w0x,w0y,wzx,wzy,vx,vy,dylaser,dzlaser,mu21,h,z0,vz,zrx,lam,cg,Ecut))
                
                prob=(1+XY['y'][2,:])/2              # Prob=(1+y[2])/2
                # probtot=probtot+prob(end)      # summing probability molecule excited at end of trajectory over all trajectories
                probtot = probtot+prob[-1]
        #   end                               # end of trajectory (skimmer 2)
      # end                                 # end of trajectory (skimmer 1)
#   end                                   # end of ntraj loop
    probtot=probtot/(icount)             # averaging probability over all trajectories
    pop.append(probtot)                   # making probability an array to plot later
    power.append(1000*P)                    # making power an array to plot later
# end                                    # end of power loop

plt.plot(power,pop)                 # plots fluence curve
plt.show()
plt.close()

if save:
    pop = np.array(pop)
    power = np.array(power)
    savedata = np.column_stack((power,pop))
    np.savetxt(savefolder+'f_'+str(fx)+'_z_'+str(z0)+'.txt',savedata, header='power, population, nmax='+str(nmax))



