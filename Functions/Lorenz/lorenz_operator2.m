function U = lorenz_operator2( tet, n )
sig_1 = [ 0 1;
          1 0 ];
sig_2 = [ 0 -1i;
          1i 0 ];
sig_3 = [ 1 0;
          0 -1 ];
U = expm(-tet / 2 * (n(1) * sig_1 + n(2) * sig_2 + n(3) * sig_3));

end

