
clear
clc
close
addpath('../Protocols');
addpath('Kraus');
addpath('Functions');
addpath('Functions/Density Matrix');
addpath('Functions/Lorenz');

fig = figure('Name','Bloch Sphere','pos',[800 200 700 700]);
figure(fig);

filename_mat = {
                'X_Tetrahedron'
%                 'X_Cube'
%                 'X_Octahedron'
%                 'X_Dodecahedron'
%                 'X_Icosahedron'
%                 'X_Fullerene'
%                 'X_Fullerene_dual'
                };
N_L = length(filename_mat);
X_matrix = {};
for k=1:N_L
    S = load(char(filename_mat(k)));
    X_matrix{k} = S.X;
end


%% Задание основных пораметров
% Радиус Блоха
global r1
% r1 = 0.90;
r1 = 1
% Число испытаний
n = 1e6;
T = 0; % Время измерений
T1 = 1; % Время амлитудной релаксации
T2 = 2; % Время фазовой релаксации
flag_a_r = 0;
flag_p_r = 0;

% Число точек на сфере Блоха
N = 80;
% Настройка шрифтов
% set(0,'defaultAxesFontSize',11)
% Отображение графиков на весь экран

%% Сканирование сферы Блоха
[X,Y,Z] = sphere(N);
type_cmap = 'jet';
cmap = colormap(type_cmap);

%% Цикл по всем протоколам
for k=1:N_L
tic
X1 = cell2mat(X_matrix(k));
dm = build_dm(0.99, pi/4, 5*pi/3);
% [R, TET, PHI] = return_r_tet_phi_by_dm(X1(1,:)'*X1(1,:));
% dm = build_dm(0.99, TET, PHI);
% dm = build_orthogonal_dm(dm);
% X1 = new_protocol_by_lorenz_transformation(X1, dm);

L_X_Cube = zeros(N + 1, N + 1);

[m, s] = size(X1);
t = time_for_protocol(n, s, m, 'uniform');

for ii = 1:size(X,1)
    for jj = 1:size(X,2)
        [phi1,tet1] = cart2sph(X(ii,jj),Y(ii,jj),Z(ii,jj));
        tet1 = pi/2 - tet1;
        dm_tmp = build_dm(r1, tet1, phi1);
%         L_X_Cube(ii,jj) =  n*real(dF_mixed_with_noise(X1, dm_tmp, t, T, T1, T2, flag_a_r, flag_p_r, n));
        L_X_Cube(ii,jj) =  real(efficiency(X1, dm_tmp, t, T, T1, T2, flag_a_r, flag_p_r, n));
    end
end

%% Закрашивание сфер Блоха L_X_Cube
% C_Cube = color
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

%% Построение сферы Блоха L_X_Cube

% surf(X,Y,Z,C_Cube,'LineStyle','none');
hold on
surf(X*r1,Y*r1,Z*r1,C_Cube,'EdgeAlpha',0.05);
hold off
colorbar;
colormap(type_cmap)
caxis([minC_Cube maxC_Cube]);
xlabel('$$x$$', 'Interpreter', 'latex');
ylabel('$$y$$', 'Interpreter', 'latex');
zlabel('$$z$$', 'Interpreter', 'latex');

title(sprintf('%s: r = %.3f, Lmin = %.3f, Lmax = %.3f, T = %.2f, T1 = %.1f', char(filename_mat(k)), r1, minC_Cube, maxC_Cube, T, T1),'Interpreter', 'none');

% Результаты
fprintf('%s\n', char(filename_mat(k)));
fprintf('Lmin = %.3f, Lmax = %.3f\n', minC_Cube, maxC_Cube);
toc
fprintf('\n');

hold on
plot_protocol_on_sphere(X1*sqrt(r1), 100, 'r');
hold off
view(3)
end
