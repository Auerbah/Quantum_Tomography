function dm = create_dm(phi, tet, p_1)
    dm_0 = build_dm(1, tet, phi);
    dm_1 = build_orthogonal_dm(dm_0);
    dm = p_1  * dm_1 + (1-p_1) * dm_0;
end

