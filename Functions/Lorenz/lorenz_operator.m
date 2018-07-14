function U = lorenz_operator( tet, tet_n, fi_n )
%% Оператор поворота
% tet - угол поворота
% tet_n - зенитный угол 
% fi_n - азимутальный угол
n1 = sin(tet_n) * cos(fi_n);
n2 = sin(tet_n) * sin(fi_n);
n3 = cos(tet_n);
sig_1 = [ 0 1;
          1 0 ];
sig_2 = [ 0 -1i;
          1i 0 ];
sig_3 = [ 1 0;
          0 -1 ];
U = expm(-tet / 2 * (n1 * sig_1 + n2 * sig_2 + n3 * sig_3));

end

