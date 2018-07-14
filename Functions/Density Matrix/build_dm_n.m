function dm = build_dm_n(n)

    I = eye(2);
    Gx = [0 1; 1 0];
    Gy = [0 -1i; 1i 0];
    Gz = [1 0; 0 -1];

    dm = I/2 + (n(1)*Gx + n(2)*Gy + n(3) *Gz)/2;

end

