%Parameter Definition
clear; close all; clc;

p=Params();
dt=0.001;  %time step [s]
Tf=4; 
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

%Heart Path

cx=0.15; cy=0.15; %Centre[m]
w=2*pi/Tf;
tau=w*t;
s=0.003; %scale
xr=cx+s*(12*sin(tau)-4*sin(3*tau));
xrd=s*w*(12*cos(tau)-12*cos(3*tau));
xrdd=s*w^2*(-12*sin(tau)+36*sin(3*tau));
xrddd=s*w^3*(-12*cos(tau)+108*cos(3*tau));
yr=cy+s*(13*cos(tau)-5*cos(2*tau)-2*cos(3*tau)-cos(4*tau)+2.5);
yrd=s*w*(-13*sin(tau)+10*sin(2*tau)+6*sin(3*tau)+4*sin(4*tau));
yrdd=s*w^2*(-13*cos(tau)+20*cos(2*tau)+18*cos(3*tau)+16*cos(4*tau));
yrddd=s*w^3*(13*sin(tau)-40*sin(2*tau)-54*sin(3*tau)-64*sin(4*tau));

%decoupled u's
ux=axisVoltage(xrd,xrdd,xrddd,p,p.bx);
uy=axisVoltage(yrd,yrdd,yrddd,p,p.by);
[ul,ur]=voltagesLR(ux,uy);

%initial states
x0=initialState(xr(1),xrd(1),xrdd(1),p,p.bx); %function below
y0=initialState(yr(1),yrd(1),yrdd(1),p,p.by);

%simulation (infinity)
zx=lsim(sys_ax,ux,t,x0); %(ss, input, time)
zy=lsim(sys_ay,uy,t,y0);

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
p.mur= 0.05; %roll frict [Ns/m]
p.mulx=0.05; %lin bear frict x [Ns/m]
p.muly=0.05; %lin bear frict y [Ns/m]
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