clear; close all; clc;
time = 0:20/1999:20;

C = false;
T_FP = 20*randi([int16(1999*15/20) int16(1999*22/20)])/1999;
T_fr = rand();
F_r = true;
T_start = 20*randi([50 100])/1999;
q_printer = 1;

T_C = time(1,1500);


for T = time

    if T == T_start
        fprintf("Print enabled\n")
        q_printer = 2;
       
    end
    
    if q_printer == 2
        q_printer = 3;
        fprintf("Printing\n")

    end
     
    if C == true && T == T_C
        q_printer = 4
        fprintf("Filament colour needs switching out\n")
        pause(3)
        fprintf("Filament colour switched out\n")
        q_printer = 2 
    end 

    if T == T_FP
    fprintf("Filament needs replacing\n")
    pause(10)
    fprintf("Filament is replaced\n")
    q_printer = 2;
    
    end
    
    if T == 20
       fprintf("Print is done")
       q_printer = 1;
    end

end


J = 0.15;



fprintf("Moving to (%4.2f, 0)", J)


