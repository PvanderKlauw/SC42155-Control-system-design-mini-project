%Parameter Definition
clear; close all; clc;

p=Params();
Tf=2; %final time [s]
dt=0.001;  %time step [s]
t=(0:dt:Tf)'; %0-Tf with dt time steps

%decoupled model
Aax=[0 1 0;
     0 0 p.rm*p.Km/p.Jeq;
     0 -p.Km/(p.rm*p.Lm) -p.Rm/p.Lm];
Bax=[0; 0; 1/p.Lm];
Cax=eye(3); %outputs 3x3 identity matrix
Dax=zeros(3,1); %outputs 3x1 column of zeros
sys_ax=ss(Aax,Bax,Cax,Dax); %makes the statespace

%Infinity sign Path
w=2*pi*0.5;  %ang freq[rad/s]
cx=0.15; cy=0.15; %Centre[m]
a=0.05; b=0.05; %Amplitudes change for elipse vs square
k=0.12; %to dictate how square/eliptical the shape is
xr=cx+a*(cos(w*t)-k*cos(3*w*t));
xrd=a*w*(-sin(w*t)+3*k*sin(3*w*t));
xrdd=a*w^2*(-cos(w*t)+9*k*cos(3*w*t));
xrddd=a*w^3*(sin(w*t)-27*k*sin(3*w*t));
yr=cy+b*(sin(w*t)+k*sin(3*w*t));
yrd=b*w*(cos(w*t)+3*k*cos(3*w*t));
yrdd=b*w^2*(-sin(w*t)-9*k*sin(3*w*t));
yrddd=b*w^3*(-cos(w*t)-27*k*cos(3*w*t));

%decoupled u's
ux=axisVoltage(xrd,xrdd,xrddd,p);
uy=axisVoltage(yrd,yrdd,yrddd,p);
[ul,ur]=voltagesLR(ux,uy);

%initial states
x0=initialState(xr(1),xrd(1),xrdd(1),p); %function below
y0=initialState(yr(1),yrd(1),yrdd(1),p);

%simulation (infinity)
zx=lsim(sys_ax,ux,t,x0); %(ss, input, time)
zy=lsim(sys_ax,uy,t,y0);

%plotting
figure;
subplot(2,1,1);
plot(t,ul,t,ur); ylabel('u [V]'); legend('u_l','u_r'); grid on;
subplot(2,1,2);
plot(t,zx(:,1),t,zy(:,1)); ylabel('Position [m]'); xlabel('t [s]'); legend('x(t)','y(t)'); grid on;

figure;
plot(zx(:,1),zy(:,1)); axis equal; grid on;
xlabel('x [m]'); ylabel('y [m]');

function p = Params()

p.rm=12e-3; %Motor pulley radius [m]
p.Lm=11e-3; %Motor inductance [H]
p.Rm=80;    %Motor resistance [Ohms]
p.Jm=9e-6;  %Motor inertia [kg/m^2]
p.Km=0.08;  %Motor constant [Nm/A]
p.m=0.05;   %mass of nozzle [kg]
p.dx=0.30;  %change of x [m]
p.dy=0.30;  %change in y [m]
p.Jeq=p.Jm+p.m*p.rm^2/2; %for subbing in Jeq

end

function u = axisVoltage(xd,xdd,xddd,p)
    u=p.Lm*p.Jeq/(p.rm*p.Km)*xddd ...
    +p.Rm*p.Jeq/(p.rm*p.Km)*xdd...
    +p.Km/p.rm*xd;
end

function [ul,ur] = voltagesLR(ux,uy)
    ul=ux+uy;
    ur=ux-uy;
end

function z0=initialState(x,xd,xdd,p)
    i0=p.Jeq*xdd/(p.rm*p.Km);
    z0=[x; xd; i0];
end