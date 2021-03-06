clc
clear
close
addpath('../Kraus')
addpath('../../Protocols');
addpath('..')

sig_1 = [ 0 1;
          1 0 ];
sig_2 = [ 0 -1i;
          1i 0 ];
sig_3 = [ 1 0;
          0 -1 ];

%% �������
fig = figure('Name','Bloch Sphere','pos',[700 200 600 600]);
figure(fig);
view(3)
xlabel('x'); ylabel('y'); zlabel('z');
xlim([-1 1]); ylim([-1 1]); zlim([-1 1]);

%% ��������� ���������
x = X_Cube;
% x = x*rotation_operator(pi/4,pi/4,pi/4);
N = length(x);
radius = zeros(1,N); tet = zeros(1,N); phi = zeros(1,N);

r = zeros(N,3);
dm = zeros(2, 2, N);

for i=1:N
    dm(:,:,i) = x(i,:)'*x(i,:);
    [radius(i), tet(i), phi(i)] = return_r_tet_phi_by_po_matrix(dm(:,:,i));
    
    r(i,:) = [radius(i)*cos(phi(i))*sin(tet(i)),...
              radius(i)*sin(phi(i))*sin(tet(i)),...
              radius(i)*cos(tet(i))];
end



%% ������������
type_cmap = 'white';
cmap = colormap(type_cmap);
hold on
[X,Y,Z] = sphere(30);
surf(X, Y, Z, 'EdgeAlpha', 0.1, 'FaceAlpha', 0.1);

for i=1:N
    scatter3(r(i,1), r(i,2), r(i,3),100, [0 0 0], 'filled')
end

%% Lorenz transformation
psi_0 = [1; 0];
psi_1 = [0; 1];

dm_0 = 0.999*psi_0*psi_0'+(1-0.999)*psi_1*psi_1';

P0 = trace(dm_0);
P1 = trace(dm_0 * sig_1);
P2 = trace(dm_0 * sig_2);
P3 = trace(dm_0 * sig_3);
velocity = [P1 P2 P3]/P0;
n = velocity/norm(velocity);
norm(velocity)
eta = atanh(norm(velocity));
[radius_0, tet_0, phi_0] = return_r_tet_phi_by_po_matrix(psi_0*psi_0');
r_0 = [radius_0*cos(phi_0)*sin(tet_0),...
       radius_0*sin(phi_0)*sin(tet_0),...
       radius_0*cos(tet_0)];  
scatter3(r_0(1), r_0(2), r_0(3), 100, 'b');
radius_1 = radius_0;
tet_1 = pi - tet_0;
phi_1 = phi_0 + pi;
r_1 = [radius_1*cos(phi_1)*sin(tet_1),...
       radius_1*sin(phi_1)*sin(tet_1),...
       radius_1*cos(tet_1)];  
scatter3(r_1(1), r_1(2), r_1(3), 100, 'b', 'filled');
L = lorenz_operator2( eta, n);
pause(1)
K = 1;

for k=1:K
for i=1:N
dm;
dm(:,:,i) = L*dm(:,:,i)*L';
dm(:,:,i) = dm(:,:,i)/trace(dm(:,:,i));
dm;
[radius(i), tet(i), phi(i)] = return_r_tet_phi_by_po_matrix(dm(:,:,i));

r(i,:) = [radius(i)*cos(phi(i))*sin(tet(i)),...
          radius(i)*sin(phi(i))*sin(tet(i)),...
          radius(i)*cos(tet(i))];

scatter3(r(i,1), r(i,2), r(i,3), 100, 'r', 'filled');
end

pause(0.2)
end
hold off
