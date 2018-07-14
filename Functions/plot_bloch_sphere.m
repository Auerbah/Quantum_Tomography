function plot_bloch_sphere(name_of_protocol, L_X_Cube, X, Y, Z, r1, type_cmap)
cmap = colormap(type_cmap);
maxC_Cube = max(max(L_X_Cube));
minC_Cube = min(min(L_X_Cube));
L_X_Cube = round((L_X_Cube - minC_Cube) / (maxC_Cube - minC_Cube) * size(cmap,1));
L_X_Cube(L_X_Cube==0) = 1;
C_Cube = zeros([size(L_X_Cube) 3]);
for i = 1:size(L_X_Cube,1)
    for j = 1:size(L_X_Cube,2)
        color = cmap(L_X_Cube(i,j),:);
        C_Cube(i,j,1) = color(1);
        C_Cube(i,j,2) = color(2);
        C_Cube(i,j,3) = color(3);
    end
end

surf(X*r1,Y*r1,Z*r1,C_Cube,'EdgeAlpha',0.05);
colorbar;
colormap(type_cmap)
caxis([minC_Cube maxC_Cube]);
xlabel('x');
ylabel('y');
zlabel('z');
title(sprintf('%s protocol:\nr = %.3f, L_{min} = %.3f, L_{max} = %.3f', name_of_protocol, r1, minC_Cube, maxC_Cube),...
    'Interpreter', 'TEX',...
    'FontSize', 16);
% title(sprintf('%s protocol:\nr = %.3f, eff_{min} = %.3f, eff_{max} = %.3f', name_of_protocol, r1, minC_Cube, maxC_Cube),...
%     'Interpreter', 'TEX',...
%     'FontSize', 16);

end

