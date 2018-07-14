function dm_new = tomography_for_weak_measurements(dm, dm_0, t_0)
    dm_1 = build_orthogonal_dm(dm_0);
    lambda_0 = fidelity(dm, dm_1);
    lambda_1 = fidelity(dm, dm_0);
    K_0 = generate_K(lambda_0, t_0, 0);
    K_1 = generate_K(lambda_1, t_0, 0);
  
    if K_0 == 0
%         K_0
        K_0 = 1;
    end
    if K_1 == 0
%         K_1
        K_1 = 1;
    end
    p_0 = K_0/K_1;
    
    dm_new = p_0*dm_1+(1-p_0)*dm_0;

end

