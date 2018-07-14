function psi = psi_phi_tet(phi,tet)

% psi = [cos(tet / 2).* exp(-1i * phi/2); sin(tet / 2) .* exp(1i * phi/2)];
psi = [cos(tet / 2); sin(tet / 2) .* exp(1i * phi)];
% psi = psi .* exp(1i * phi/2);
