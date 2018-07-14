function [phi, tet] = return_phi_tet_2(psi)

phi = phase(psi(2))-phase(psi(1));
if phi < 0 
    phi = phi + 2 * pi;
end
    
tet = 2*acos(abs(psi(1)));
