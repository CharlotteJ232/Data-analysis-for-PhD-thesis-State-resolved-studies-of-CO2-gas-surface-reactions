function dydx = bloch_diffequ(t,y,E0,w0x,w0y,wzx,wzy,vx,vy,dylaser,dzlaser,mu21,h,z0,vz,zrx,lambda,cg,Ecut)

%%%%%%%%%%%the Optical Bloch equation
dydx = zeros(3,1);

if (E0*(exp(-vx*t*vx*t/(wzx*wzx)-((dylaser+vy*t)*(dylaser+vy*t)/(wzy*wzy)))) > Ecut);
dydx(1)=-2*pi*(vz/lambda-vx*vx*t/(lambda*(dzlaser+vz*t+z0+(zrx*zrx)/(dzlaser+vz*t+z0))))*y(2);
dydx(2)= 2*pi*((vz/lambda-vx*vx*t/(lambda*(dzlaser+vz*t+z0+(zrx*zrx)/(dzlaser+vz*t+z0))))*y(1)+(cg*mu21*E0/h)*exp(-vx*t*vx*t/((w0x*sqrt(1+(dzlaser+vz*t+z0)*(dzlaser+vz*t+z0)/(zrx*zrx)))*(w0x*sqrt(1+(dzlaser+vz*t+z0)*(dzlaser+vz*t+z0)/(zrx*zrx))))-(dylaser+vy*t)*(dylaser+vy*t)/(wzy*wzy))*sqrt((w0x*w0y)/((w0x*sqrt(1+(dzlaser+vz*t+z0)*(dzlaser+vz*t+z0)/(zrx*zrx)))*wzy))*y(3));
dydx(3)=-2*pi*((cg*mu21*E0/h)*exp(-vx*t*vx*t/((w0x*sqrt(1+(dzlaser+vz*t+z0)*(dzlaser+vz*t+z0)/(zrx*zrx)))*(w0x*sqrt(1+(dzlaser+vz*t+z0)*(dzlaser+vz*t+z0)/(zrx*zrx))))-((dylaser+vy*t)*(dylaser+vy*t)/(wzy*wzy)))*sqrt((w0x*w0y)/((w0x*sqrt(1+(dzlaser+vz*t+z0)*(dzlaser+vz*t+z0)/(zrx*zrx)))*wzy)))*y(2);
end

% x=vx*t
% y=vy*t+dylaser
% z=dzlaser+vz*t+z0

% Rz=(dzlaser+vz*t+z0+(zrx*zrx)/(dzlaser+vz*t+z0))=z+zrx*zrx/z
% wzx=((w0x*sqrt(1+(dzlaser+vz*t+z0)*(dzlaser+vz*t+z0)/(zrx*zrx))=w0x*sqrt(1+z*z/(zrx/zrx))

% Frequencies in rad/s, and Bloch equations are:
% dy(1)/dt=-(doppler-sweep)*y(2)
% dy(2)/dt=(doppler-sweep)*y(1)+omega(z)*y(3)
% dy(3)/dt=-omega(z)*y(2)

% doppler=vz/lambda
% sweep=vx*vx*t/(Rz*lambda)
% omega(z)=(mu*E/h)*sqrt(wx_0*wy_0/(wx_z*wy_z)
% E=E_0*exp(-x*x/(wx_z*wx_z)-y*y/(wy_z*wy_z))