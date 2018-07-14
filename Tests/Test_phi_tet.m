clc
clear
close all

addpath('../Kraus')
addpath('..')

global d
d = 1e-10

psi0 = [1; 0]; % {0, 0, 1}
[phi, tet] = return_phi_tet_2(psi0);
equals_phi_tet(phi, 0, tet, 0);
psi = psi_phi_tet(phi,tet);
equals_psi(psi, psi0);

[r, tet1, phi1] = return_r_tet_phi_by_po_matrix(psi0*psi0');
equals_phi_tet(phi1, 0, tet1, 0);

psi1 = [0; 1]; %{0, 0, -1}
[phi, tet] = return_phi_tet_2(psi1);
equals_phi_tet(phi, 0, tet, pi)
psi = psi_phi_tet(phi,tet);
equals_psi(psi, psi1);

[r, tet1, phi1] = return_r_tet_phi_by_po_matrix(psi1*psi1');
equals_phi_tet(phi1, 0, tet1, pi)

phi0 = 1 / sqrt(2) * [1; 1]; %{1, 0 ,0}
[phi, tet] = return_phi_tet_2(phi0);
equals_phi_tet(phi, 0, tet, pi/2)
psi = psi_phi_tet(phi,tet);
equals_psi(psi, phi0);

[r, tet1, phi1] = return_r_tet_phi_by_po_matrix(phi0*phi0');
equals_phi_tet(phi1, 0, tet1, pi/2)

phi1 = 1 / sqrt(2) * [1; -1];
[phi, tet] = return_phi_tet_2(phi1);
equals_phi_tet(phi, pi, tet, pi/2)
psi = psi_phi_tet(phi,tet);
equals_psi(psi, phi1);

[r, tet1, phi1] = return_r_tet_phi_by_po_matrix(phi1*phi1');
equals_phi_tet(phi1, pi, tet1, pi/2)

chi0 = 1 / sqrt(2) * [1; 1i];
[phi, tet] = return_phi_tet_2(chi0);
equals_phi_tet(phi, pi/2, tet, pi/2)
psi = psi_phi_tet(phi,tet);
equals_psi(psi, chi0);

[r, tet1, phi1] = return_r_tet_phi_by_po_matrix(chi0*chi0');
equals_phi_tet(phi1, pi/2, tet1, pi/2)

chi1 = 1 / sqrt(2) * [1; -1i];
[phi, tet] = return_phi_tet_2(chi1);
equals_phi_tet(phi, 3*pi/2, tet, pi/2)
psi = psi_phi_tet(phi,tet);
equals_psi(psi, chi1);

[r, tet1, phi1] = return_r_tet_phi_by_po_matrix(chi1*chi1');
equals_phi_tet(phi1, 3*pi/2, tet1, pi/2)

function equals_phi_tet(phi0, phi1, tet0, tet1)
    global d
    if abs(phi0-phi1)>d || abs(tet0-tet1)>d
        fprintf("phi %f != %f\nor\n%f != %f \n",phi0,phi1,tet0,tet1);
        error('not equal') 
    end
end

function equals_psi(psi0, psi1)
    global d
    res = psi0 - psi1;
    for i=1:length(psi0)
    if abs(psi0(i)-psi1(i))>d
        psi0
        psi1
        error('not equal')
    end
    end
end