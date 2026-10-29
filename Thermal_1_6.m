%Ellipse
%Parameter Definition
clear; close all; clc;

p=Params();
dt=0.01;  %time step [s]
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
a=0.05; b=0.05; %Amplitudes change for elipse vs square
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
zy=lsim(sys_ay,uy,t,y0);

%start of 1.6
vel = sqrt(zx(:,2).^2 + zy(:,2).^2); %vel mag
s_path = cumtrapz(t, vel); %cumulative distance along path
L_total = s_path(end); %total length of printed ellipse

%Thermal Discretization
dx = 0.001; %spatial resolution [m]
N_nodes = ceil(L_total / dx) + 1;%total spatial nodes
T = NaN(N_nodes, 1); %Initialize temperature array (NaN = unprinted)

%Thermal diffusivity
alpha = p.k / (p.rho * p.cp);
C1 = alpha * dt / dx^2;
C2 = (p.hb / (p.rho * p.cp * p.a)) * dt;
C3 = (3 * p.hair / (p.rho * p.cp * p.a)) * dt;

%Target steady-state uniform temperature
T_eq = (p.hb * p.Tb + 3 * p.hair * p.Te) / (p.hb + 3 * p.hair);
uniform_tolerance = 0.1; % [C]

Phi=zeros(N_nodes+1, 1);
f_T=zeros(N_nodes, 1);

%Active Printing 
for n = 1:(length(t)-1)
    current_length = s_path(n); %nodes currently active based on distance traveled
    active_idx = floor(current_length / dx) + 1;
    T(active_idx) = p.Tn; %BCs
    T_new = T;
    f_T = C2*(p.Tb - T) - C3*(T - p.Te);
    
    if active_idx > 1 
        Phi(2:active_idx) = (T(1:active_idx-1) - T(2:active_idx)) * C1;
        Phi(1) = -Phi(2); %boundary flow
        T_new(1) = T(1) + Phi(1) - Phi(2) + f_T(1); %(Phi_in-Phi_out+sources)
    end
    
    if active_idx > 2
        %interior nodes (Phi_in - Phi_out + sources)
        i = 2:(active_idx-1);
        T_new(i) = T(i) + Phi(i) - Phi(i+1) + f_T(i);
    end
    
    T = T_new;
end

T(isnan(T))=p.Tn;
% Post-Print Cooling
t_cool = t(end);
cooling = true;

while cooling
    t_cool = t_cool + dt;
    Phi(2:N_nodes) = (T(1:N_nodes-1) - T(2:N_nodes)) * C1;
    
    %Boundary flows
    Phi_periodic = (T(N_nodes) - T(1)) * C1;
    Phi(1) = Phi_periodic;
    Phi(N_nodes+1) = Phi_periodic;
    
    %distributed sources
    f_T = C2*(p.Tb - T) - C3*(T - p.Te);
    T = T + Phi(1:N_nodes) - Phi(2:N_nodes+1) + f_T;
    max_deviation = max(abs(T - T_eq));
    if max_deviation < uniform_tolerance
        cooling = false;
    end
end
fprintf('Total required cooling time (Printing + Post-Cooling): %.2f seconds\n', t_cool);
fprintf('Uniform equilibrium temperature reached: %.2f C\n', T_eq);

%plotting
figure;
subplot(2,1,1);
plot(t,ul,t,ur); ylabel('u [V]'); legend('u_l','u_r'); grid on;
subplot(2,1,2);
plot(t,zx(:,1),t,zy(:,1)); ylabel('Position [m]'); xlabel('t [s]'); legend('x(t)','y(t)'); grid on;

figure;
plot(zx(:,1),zy(:,1)); axis equal; grid on;
xlabel('x [m]'); ylabel('y [m]');

figure;
x_axis = linspace(0, L_total, N_nodes);
plot(x_axis, T, 'LineWidth', 1.5);
yline(T_eq, 'r--', 'Steady State T_{eq}');
xlabel('Position along Filament Path [m]');
ylabel('Temperature [C]');
title('Final Filament Temperature Distribution');
grid on;

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

%Thermal Params
p.rho=1250; %Density (PLA) [kg/m^3]
p.cp=1800; %Specfic heat (PLA) [J/(kgK)]
p.k = 0.13; %Thermal conductivity of PLA [W/(mK)]
p.a = 0.4e-3; %Square cross section side length [m]
p.Tn = 210; %Nozzle extrusion temperature [C]
p.Tb = 60; %Bed temperature [C]
p.Te = 20; %Ambient air temperature [C]
p.hair = 15; %Natural convection coefficient [W/(m^2K)]
p.hb = 100; %Contact conductance [W/(m^2K)]
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