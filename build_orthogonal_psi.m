function psi = build_orthogonal_psi(psi)
    dm = psi*psi';
    dm = build_orthogonal_dm(dm)
    [r, tet, phi] = return_r_tet_phi_by_dm(dm);
    
    psi = psi_phi_tet(phi,tet)
end