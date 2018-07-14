clc
clear
close all

hold on 

[X,Y,Z] = sphere(30);


s1 = surf(X, Y, Z, 'EdgeAlpha', 0.1, 'FaceAlpha', 0.1);
cmap1 = colormap('winter');

s2 = surf(0.5*X, 0.5*Y, 0.5*Z, 'EdgeAlpha', 0.1, 'FaceAlpha', 1);
hold off

view(3)