function [r, tet, phi] = return_r_tet_phi_by_dm(po)
% Функция возвращает r, tet, phi для заданной матрицы плотности
% Для этого используется представление:
% po = I/2 + (nx*Gx + ny*Gy + nz *Gz)/2
% или
% po = 1/2 * [1+nz nx-1i*ny;
%             nx+1i*ny 1-nz];
% где nx = r*cos(phi)*sin(tet);
%     ny = r*sin(phi)*sin(tet);
%     nz = r*cos(tet);

    nx = real(po(2,1)+po(1,2));
    ny = real(1i*(po(1,2)-po(2,1)));
    nz = real(po(1,1)-po(2,2));
    
    [phi,tet,r] = cart2sph(nx,ny,nz);
    
    tet = pi/2 - tet;
end