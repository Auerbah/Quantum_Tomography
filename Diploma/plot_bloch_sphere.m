clc
clear
close
addpath('../Functions/Density Matrix')

%% Графика
fig = figure('Name','Bloch Sphere','pos',[700 200 800 800]);
figure(fig);
view(3)
% Set the 'visible' property 'off'
set(gca,'XTick',[])
set(gca,'YTick',[])
set(gca,'ZTick',[])
%% Начальное состояние

psi_0 = [1;0];
psi_1 = [0;1];

dm_0 = psi_0*psi_0';
dm_1 = psi_1*psi_1';
dm_2 = 0.7*dm_0 + 0.3*dm_1;

dm = zeros(2,2,2);
dm(:,:,1) = dm_0;
dm(:,:,2) = dm_2;

N = length(dm);

radius = zeros(1,N);
tet = zeros(1,N);
phi = zeros(1,N);
r = zeros(N,3);
for i=1:N
    [radius(i), tet(i), phi(i)] = return_r_tet_phi_by_dm(dm(:,:,i));
    
    r(i,:) = [radius(i)*cos(phi(i))*sin(tet(i)),...
              radius(i)*sin(phi(i))*sin(tet(i)),...
              radius(i)*cos(tet(i))];
end

rad = radius(2);

%% Визуализация
hold on
[X,Y,Z] = sphere(30);
colormap(fig, winter);

surf1 = surf(X, Y, Z, 'EdgeAlpha', 0, 'FaceAlpha', 0.1);

[X,Y,Z] = sphere(90);
surf2 = surf(rad*X, rad*Y, rad*Z, 'EdgeAlpha', 0, 'FaceAlpha', 0.5);

scatter3(r(1,1), r(1,2), r(1,3), 50, [0 0 0], 'filled')
scatter3(r(2,1), r(2,2), r(2,3), 100, [0 0 0], 'filled')



hold off