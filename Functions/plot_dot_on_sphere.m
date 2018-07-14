function plot_dot_on_sphere(dm, size, color)
[r, tet, phi] = return_r_tet_phi_by_dm(dm);
coordinates = [r*cos(phi)*sin(tet),...
               r*sin(phi)*sin(tet),...
               r*cos(tet)];
scatter3(coordinates(:,1), coordinates(:,2), coordinates(:,3), size, color, 'filled')
end

