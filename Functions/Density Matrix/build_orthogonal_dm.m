function dm = build_orthogonal_dm(dm)
    [r, tet, phi] = return_r_tet_phi_by_dm(dm);
    dm = build_po_matrix(r, pi-tet, pi+phi);
end

