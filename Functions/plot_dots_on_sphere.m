function plot_dots_on_sphere(dm, size, color)
N = length(dm);
for i = 1:N
[r, tet, phi] = return_r_tet_phi_by_po_matrix(dm(:,:,i));
coordinates = [r*cos(phi)*sin(tet),...
               r*sin(phi)*sin(tet),...
               r*cos(tet)];
scatter3(coordinates(:,1), coordinates(:,2), coordinates(:,3), size, color, 'filled')
end
end



