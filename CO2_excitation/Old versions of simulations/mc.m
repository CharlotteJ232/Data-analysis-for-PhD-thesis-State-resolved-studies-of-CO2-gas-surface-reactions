clear
close
clc
% To simulate power fluence curves with gaussian velocities etc

%
% Laser beam defined as travelling along z axis, molecular beam defined
% as travelling along x axis.
%

%% Expt input variables
lambda=3.30e-6; % wavelength (m)
fx=0.254;       % focus length of lens in xz plane (m)
rx=2e-03;       % radius of laser in x (m)
ry=2e-03;       % radius of laser in y direction (m)
vel=2056;       % velocity of molecular beam (m/s)
fwhmv=0.1;      % full width half maximum of velocity distribution (as fraction of v)
ang=0.35;       % angular spread of beam out of nozzle (degrees)

z0=0.106;   % position of molecular beam from focussed waist (m)
A21=25.35;  % Einstein coefficient (Hz)
cg=sqrt(1); % transition probability (CG coefficient)

pcut=5e-05; % minimum laser power for integration (W)
pmin=0;     % minimum laser power for fluence curve (W)
pmax=1;     % maximum laser power for fluence curve (W)
pstep=0.02; % laser power step for fluence curve (W)

nmax=500;   % number of trajectories

noztoskim=60e-03;   % distance between nozzle and 1st skimmer (m)
noztoskim2=160e-03; % distance between nozzle and 2nd skimmer (m)
noztolaser=380e-03; % distance between nozzle and laser (m)

znozzle=0; % offset of nozzle to skimmers along laser axis (m)
ynozzle=0; % offset of nozzle to skimmers perpendicular to laser axis (m)

skimdiam=1e-03;          % 1st skimmer diameter (m)
skimdiam=skimdiam/2;     % 1st skimmer radius (m)
skimdiam2=3e-03;         % 2nd skimmer diameter (m)
skimdiam2=skimdiam2/2;   % 2nd skimmer radius (m) 
%% constants
eps0=8.85e-12;	% vacuum permitivity (F/m)
h=6.626e-34;    % Plancks constant (Js)
c= 3e08;        % speed of light (m/s)
ip=0;           % counter for powers

%% calculating values that are the same for each trajectory
w0x=lambda*fx/(pi*rx);   % x beam waist at focal point of laser (m)
w0y=ry;                  % laser beam not focussed in y direction
zrx=pi*w0x*w0x/lambda;   % Rayleigh range before focussing (m)
A=pi*w0x*w0y;            % area of focussed laser beam at waist (m^2)
wb=2*pi*c/lambda;        % transition frequency
mu21=sqrt((3*eps0*h*c*c*c*A21)/(2*(wb*wb*wb))); % dipole moment (Cm)

%% Looping over trajectories
 for P=pmin:pstep:pmax      % looping over laser power
  ip=ip+1                   % Counter for powers used in a later loop
  I=2*P/A;                  % maximum laser intensity at waist (W/m^2)
  E0=sqrt((2*I)/(eps0*c));  % electic field at waist (V/m)
  iskim=0;                  % zeroing counters for molecules that don't make it through the skimmers
  icount=0;                 % zeroing counters for molecules that do make it through the skimmer
  probtot=0;                % zeroing counters for probability
 
  parfor ntraj=1:1:nmax;       % looping over trajectories
    % velocity selection
    sintheta=rand(1);       % polar angle of v - need to sample from sin distribution due to sampling over sphere
    theta=asin(sintheta);
    theta=theta*ang/90;
    phi=rand(1)*2*pi;         % azimuthal angle of v
    v=randn(1)*fwhmv*vel+vel; % velocity
    vz=v*sin(theta)*cos(phi); % initial velocity in z direction
    vy=v*sin(theta)*sin(phi); % initial velocity in y direction
    vx=v*cos(theta);          % initial velocity in x direction
     
    % 1st skimmer
    tskim=noztoskim/vx;                         % time to reach 1st skimmer
    dyskim=tskim*vy-ynozzle;                    % y position wrt 1st skimmer
    dzskim=tskim*vz-znozzle;                    % z position wrt 1st skimmer
    drskim=sqrt(dyskim*dyskim+dzskim*dzskim);   % determining if trajectory goes through 1st skimmer
    
    % 2nd skimmer
    tskim2=noztoskim2/vx;                       % time to reach 2nd skimmer
    dyskim2=tskim2*vy-ynozzle;                  % y position wrt 2nd skimmer
    dzskim2=tskim2*vz-znozzle;                  % z position wrt 2nd skimmer
    drskim2=sqrt(dyskim2*dyskim2+dzskim2*dzskim2); % determining if trajectory goes through 2nd skimmer

    if (drskim > skimdiam);
      iskim=iskim+1;            % molecule doesn't go through 1st skimmer
    else
      if (drskim2 > skimdiam2);  
        iskim=iskim+1;          % molecule doesn't go through 2nd skimmer
      else
        icount=icount+1;                % molecule does go through skimmer
        tlaser=noztolaser/vx;           % time taken to reach laser centre (s)
        dylaser=tlaser*vy;              % distance of molecule perpendicular to plane of laser (m)
        dzlaser=tlaser*vz;              % displacement of molecule from z0 (m)
        z=z0+dzlaser;                   % distance of molecule from focussed waist (m)
        wzx=w0x*sqrt(1+z*z/(zrx*zrx));  % beam waist at position z (m)
        wzy=ry;                         % beam not focussed in y direction (m)
        T=wzx/vx;                       % half transit time of beam
        Icut=2*pcut/(wzx*wzy);          % intensity cut off for integration
        Ecut=sqrt((2*Icut)/(eps0*c));   % electric field cut off for integration (approx 200 works well)
        [X,Y] = ode45(@(t,y) bloch_diffequ(t,y,E0,w0x,w0y,wzx,wzy,vx,vy,dylaser,dzlaser,mu21,h,z0,vz,zrx,lambda,cg,Ecut),[-3*T,3*T],[0,0,-1]);
        prob=(1+Y(:,3))/2;              % Prob=(1+y(3))/2
        probtot=probtot+prob(end);      % summing probability molecule excited at end of trajectory over all trajectories
      end                               % end of trajectory (skimmer 2)
    end                                 % end of trajectory (skimmer 1)

  end                                   % end of ntraj loop
  probtot=probtot/(icount);             % averaging probability over all trajectories
  pop(ip)=probtot;                      % making probability an array to plot later
  power(ip)=1000*P;                     % making power an array to plot later
 end                                    % end of power loop

%% Outputting power curve
%   A=xlsread('data_R1.xlsx');  % Uncomment to read in experimental data

iii=1:1:ip;                     % looping over laser power
plot(power,pop)                 % plots fluence curve

%  xlswrite('Theory.xlsx',[power',pop'])    % Uncomment to write data to Excel
%  hold on                                  % Uncomment to compare experiment and theory
%   plot(A(:,1),A(:,2),'*')                 % Plot experimental data, colon means use all A(row,column) 
%      ylim([0 1.1])                        % set y limits
%   xlim([0 1000])                          % set x limits
