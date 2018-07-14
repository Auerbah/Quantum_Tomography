clc
clear
close
t= 1;
T1 = 1;
p =0;
T2 = 1;
    Ea0 = [[1 0]
           [0 exp(-t/(2*T1))]];
    Ea0 = sqrt(p) * Ea0;
    
    Ea1 = [[0 sqrt(1-exp(-t/T1))]
          [0 0]];
    Ea1 = sqrt(p) * Ea1;
    
    Ep0 = [[1 0]
          [0 exp(-t/(2*T2))]];
    Ep0 = sqrt(1-p) * Ep0;
    Ep1 = [[0 0]
          [0 sqrt(1-exp(-t/T2))]];
    Ep1 = sqrt(1-p) * Ep1;
        Ea0'*Ea0+Ea1'*Ea1+...
        Ep0'*Ep0+Ep1'*Ep1