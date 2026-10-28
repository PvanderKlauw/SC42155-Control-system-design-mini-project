clear; close all; clc;
F_r= false;
C = false;
T = 0:20/1999:20;
T_FP = normrnd(17,0.5);
T_fr = rand();
F_r = true;
fprintf("Print enabled")




if F_r == true
    a = 0;
    b = 1;
    c = 1;
   
    r = ones(1,length(T))*normrnd(0,1.5);
    xr = b *cos(2*pi()*T/20);
    yr = c *sin(2*pi()*T/20);
    figure
    subplot(2,1,1)
    plot(yr,xr)
    subplot(2,1,2)
    plot(T,yr,T,xr)
    fprintf("Print succesful"\n)
    T_eq = ones(1,length(T))*210;

T_noz = T_eq - exp(-T/10).*r.*cos(-pi()*T);
end    






