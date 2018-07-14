clc
clear
close
addpath('../Kraus')
addpath('..')
% psi = [1;1];
% psi = psi/norm(psi)
% po = psi*psi';
% t = 2;
% T1 = 1;
% po = E_ampl(po,t, T1)
% [radius, tet, phi] = return_r_tet_phi_by_po_matrix(po)

p = 0.5
psi0 = [1; 0];
psi1 = [0; 1];

phi0 = 1 / sqrt(2) * [1; 1];
phi1 = 1 / sqrt(2) * [1; -1];

chi0 = 1 / sqrt(2) * [1; 1i];
chi1 = 1 / sqrt(2) * [1; -1i];

po = p * psi0*psi0' + (1-p) * phi1*psi1'
po = po/trace(po)
check_po_matrix(po)