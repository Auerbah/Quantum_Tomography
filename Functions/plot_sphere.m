function plot_sphere()
    az = -28;
    el = 25;
    view(az, el);
    cmap = colormap('white');
    [X,Y,Z] = sphere(30);
    surf(X, Y, Z, 'EdgeAlpha', 0.1, 'FaceAlpha', 0.1);
%     view(3)
end

