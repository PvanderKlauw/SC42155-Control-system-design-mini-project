%Ellipse
%Parameter Definition
clear; close all; clc;

p=Params();
dt=0.001;  %time step [s]
Tf=2; 
t=(0:dt:Tf);
%decoupled model
Aax=[0 1 0;
     0 -p.bx/p.Jeq p.rm*p.Km/p.Jeq;
     0 -p.Km/(p.rm*p.Lm) -p.Rm/p.Lm];
Aay=[0 1 0;
     0 -p.by/p.Jeq p.rm*p.Km/p.Jeq;
     0 -p.Km/(p.rm*p.Lm) -p.Rm/p.Lm];
Bax=[0; 0; 1/p.Lm];
Cax=eye(3); %outputs 3x3 identity matrix
Dax=zeros(3,1); %outputs 3x1 column of zeros
sys_ax=ss(Aax,Bax,Cax,Dax); %makes the statespace
sys_ay=ss(Aay,Bax,Cax,Dax);

%Ellipse Path
cx=0.15; cy=0.15; %Centre[m]
w=2*pi*0.5;  %ang freq[rad/s]
a=0.05; b=0.03; %Amplitudes change for elipse vs square
xr=cx+a*cos(w*t);
xrd=a*w*-sin(w*t);
xrdd=a*w^2*-cos(w*t);
xrddd=a*w^3*sin(w*t);
yr=cy+b*sin(w*t);
yrd=b*w*cos(w*t);
yrdd=b*w^2*-sin(w*t);
yrddd=b*w^3*-cos(w*t);

%decoupled u's
ux=axisVoltage(xrd,xrdd,xrddd,p,p.bx);
uy=axisVoltage(yrd,yrdd,yrddd,p,p.by);
[ul,ur]=voltagesLR(ux,uy);

%initial states
x0=initialState(xr(1),xrd(1),xrdd(1),p,p.bx); %function below
y0=initialState(yr(1),yrd(1),yrdd(1),p,p.by);

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

p.rm=6e-3; %Motor pulley radius [m]
p.Lm=3e-3; %Motor inductance [H]
p.Rm=3;    %Motor resistance [Ohms]
p.Jm=6e-6;  %Motor inertia [kg/m^2]
p.Km=0.05;  %Motor constant [Nm/A]
p.m=0.15;   %mass of nozzle [kg]
p.dx=0.35;  %change of x [m]
p.dy=0.35;  %change in y [m]
p.mur= 0.5; %roll frict [Ns/m]
p.mulx=0.002; %lin bear frict x [Ns/m]
p.muly=0.001; %lin bear frict y [Ns/m]
p.Jeq=p.Jm+p.m*p.rm^2/2; %for subbing in Jeq
p.b1=p.mur*p.rm^2+p.rm^2/4*(p.mulx+p.muly);
p.b2=p.rm^2/4*(p.mulx-p.muly);
p.bx=p.b1+p.b2; %x damping
p.by=p.b1-p.b2; %y damping

end

function u = axisVoltage(xd,xdd,xddd,p,b)
    u=p.Lm*p.Jeq/(p.rm*p.Km)*xddd ...
    +(p.Lm*b+p.Rm*p.Jeq)/(p.rm*p.Km)*xdd...
    +(p.Rm*b+p.Km^2)/(p.rm*p.Km)*xd;
end

function [ul,ur] = voltagesLR(ux,uy)
    ul=ux+uy;
    ur=ux-uy;
end

function z0=initialState(x,xd,xdd,p,b)
    i0=(p.Jeq*xdd+b*xd)/(p.rm*p.Km);
    z0=[x; xd; i0];
end