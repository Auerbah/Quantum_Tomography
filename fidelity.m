function F = fidelity(po1, po2)
%Fidelity 2х матриц через trace

    F = real(trace((po2^(1/2)*po1*po2^(1/2))^(1/2)))^2;
    if F > 1 + 1e-7
       F
       error('Fidelity  > 1'); 
       
    end
end

