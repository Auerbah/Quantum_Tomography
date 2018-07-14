clc
clear
close all
addpath('../Kraus')
addpath('../../Protocols');
addpath('..')

%% Графика
fig = figure('Name','Bloch Sphere','pos',[700 200 800 800]);
figure(fig);
view(3)
xlabel('x'); ylabel('y'); zlabel('z');
xlim([-1 1]); ylim([-1 1]); zlim([-1 1]);

%% Начальное состояние
x = X_Cube;
N = length(x);

radius = zeros(1,N);
% radius = ones(1,N);
tet = zeros(1,N); phi = zeros(1,N);

r = zeros(N,3);
dm = zeros(2, 2, N);
psi = x';
for i=1:N
    dm(:,:,i) = x(i,:)'*x(i,:);
    [radius(i), tet(i), phi(i)] = return_r_tet_phi_by_po_matrix(dm(:,:,i));

%     [phi(i), tet(i)] = return_phi_tet_2(psi(:,i));
    
    r(i,:) = [radius(i)*cos(phi(i))*sin(tet(i)),...
              radius(i)*sin(phi(i))*sin(tet(i)),...
              radius(i)*cos(tet(i))];
end



%% Визуализация
type_cmap = 'white';
cmap = colormap(type_cmap);
hold on
[X,Y,Z] = sphere(30);
surf(X, Y, Z, 'EdgeAlpha', 0.1, 'FaceAlpha', 0.1);

for i=1:N
    scatter3(r(i,1), r(i,2), r(i,3), 80, 'b', 'filled')
end

%% Lorenz transformation
eta = 0.5;
L = [exp(-eta/2) 0;
     0 exp(eta/2)];
tet_n = 0;
fi_n = 0;
L = lorenz_operator( eta, tet_n, fi_n);
pause(1)
for i=6
scatter3(r(i,1), r(i,2), r(i,3), 80, 'g', 'filled')
dm(:,:,i) = L*dm(:,:,i)*L';
dm(:,:,i) = dm(:,:,i)/trace(dm(:,:,i));
phi(i)
tet(i)
[radius(i), tet(i), phi(i)] = return_r_tet_phi_by_po_matrix(dm(:,:,i));
phi(i)
tet(i)
% psi(:,i) = L*psi(:,i);
% [phi(i), tet(i)] = return_phi_tet_2(psi(:,i));

r(i,:) = [radius(i)*cos(phi(i))*sin(tet(i)),...
          radius(i)*sin(phi(i))*sin(tet(i)),...
          radius(i)*cos(tet(i))];

scatter3(r(i,1), r(i,2), r(i,3), 50, 'r', 'filled');
pause(1)
end
hold off

