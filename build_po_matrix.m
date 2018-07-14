function po = build_po_matrix(r, tet, phi)
% Функция возвращает матрицу плотности по заданным r, tet, phi

    nx = r*cos(phi)*sin(tet);
    ny = r*sin(phi)*sin(tet);
    nz = r*cos(tet);

    I = eye(2);
    Gx = [0 1; 1 0];
    Gy = [0 -1i; 1i 0];
    Gz = [1 0; 0 -1];

    po = I/2 + (nx*Gx + ny*Gy + nz *Gz)/2;

end

