function L = find_lorenz_transformation(dm)
    sig_1 = [ 0 1;
              1 0 ];
    sig_2 = [ 0 -1i;
             1i 0 ];
    sig_3 = [ 1 0;
              0 -1 ];
    P0 = trace(dm);
    P1 = trace(dm * sig_1);
    P2 = trace(dm * sig_2);
    P3 = trace(dm * sig_3);
    velocity = [P1 P2 P3]/P0;
    n = velocity/norm(velocity);
    eta = atanh(norm(velocity));
    L = lorenz_operator2( eta, n);
end

