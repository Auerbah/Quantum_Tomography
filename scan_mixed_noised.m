
clear
clc
close
addpath('../Protocols');
addpath('../Protocols/Simple');
addpath('../Protocols/Belinsky');
addpath('Kraus');
addpath('../Protocols/Thomson/X');
addpath('../Noisy Measurements');
addpath('../Noisy Measurements/Scanning Bloch Sphere');

fig = figure('Name','Bloch Sphere','pos',[800 200 700 700]);
figure(fig);

filename_mat = {
%                   'X_thomson_4'
%                   'X_thomson_5'
%                   'X_thomson_6'
%                   'X_thomson_7'
%                   'X_thomson_8'
%                   'X_thomson_9'
%                   'X_thomson_10'
%                   'X_thomson_11'
%                   'X_thomson_12'
%                   'X_thomson_13'
%                   'X_thomson_14'
%                   'X_thomson_15'
%                   'X_thomson_16'
%                   'X_thomson_17'
%                   'X_thomson_18'
%                   'X_thomson_19'
%                   'X_thomson_20'
%                   'X_thomson_21'
%                   'X_thomson_22'
%                   'X_thomson_23'
%                   'X_thomson_24'
%                   'X_thomson_25'
%                   'X_thomson_26'
%                   'X_thomson_27'
%                   'X_thomson_28'
%                   'X_thomson_29'
%                   'X_thomson_30'
%                   'X_thomson_31'
%                   'X_thomson_32'
%                   'X_thomson_33'
%                   'X_thomson_34'
%                   'X_thomson_35'
%                   'X_thomson_36'
%                 'pack__1qb_4st'
%                 'pack__1qb_5st'
%                 'pack__1qb_6st'
%                 'pack__1qb_7st'
%                 'pack__1qb_8st'
%                 'pack__1qb_9st'
%                 'pack__1qb_10st'
%                 'pack__1qb_11st'
%                 'pack__1qb_12st'
%                 'pack__1qb_24st'
%                 'pack__1qb_72st'
%                 'pack__1qb_100st'
%                 'X_Tetrahedron'
%                 'X_Cube'
                'X_Octahedron'
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
r1 = 0.5
% Число испытаний
n = 4000;
T = 0.02; % Время измерений
T1 = 1; % Время амлитудной релаксации
T2 = 2; % Время фазовой релаксации
flag_a_r = 1;
flag_p_r = 1;

% Коэффициент
koeff = 1;
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
X1 = cell2mat(X_matrix(k))*rotation_operator(pi/32,pi/32,pi/32);
L_X_Cube = zeros(N + 1, N + 1);

[m, s] = size(X1);
t = time_for_protocol(n, s, m, 'uniform');

for ii = 1:size(X,1)
    for jj = 1:size(X,2)
        [phi1,tet1] = cart2sph(X(ii,jj),Y(ii,jj),Z(ii,jj));
        tet1 = pi/2 - tet1;
        po = build_po_matrix(r1, tet1, phi1);
%         if(rank(po) == 2)
%             fprintf('rank = %d\nr = %f, tet = %f, phi = %f\n', rank(po),r1,tet1,phi1);
%             L_X_Cube(ii,jj) =  koeff * n * real(dF_mixed_no_noise_new(X1, po, t, n));
%         else
            L_X_Cube(ii,jj) =  koeff * n * real(dF_mixed_with_noise(X1, po, t, T,T2,flag_a_r,flag_p_r, T1, n, rank(po)));
%             L_X_Cube(ii,jj)
%         end
    end
end

%% Закрашивание сфер Блоха L_X_Cube
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

% new_filename_png = sprintf('../Картинки/Pure All Protocols/png/Thomson/%s', char(filename_mat(k)));
% new_filename_fig = sprintf('../Картинки/Pure All Protocols/fig/Thomson/%s', char(filename_mat(k)));
% 
% savefig(fig,new_filename_fig);
% print(fig, new_filename_png, '-dpng');

%
r = [];
for i=1:length(X1)
    [phi tet] = return_phi_tet_2(conj(X1(i,:)));
    r = [r; [cos(phi)*sin(tet), sin(phi)*sin(tet), cos(tet)]];
end
r = r*r1;
hold on
scatter3(r(:,1), r(:,2), r(:,3), 80, 'r', 'filled')
hold off
view(3)

%%
%{
figure;
type_cmap = 'white';
cmap = colormap(type_cmap);
hold on
[X,Y,Z] = sphere(30);
surf(X*r1, Y*r1, Z*r1, 'EdgeAlpha', 0.1, 'FaceAlpha', 0.1);
scatter3(r(:,1), r(:,2), r(:,3), 80, 'r', 'filled')
eps = 2;
for K=1:length(r)
n = find_nearest_k(K,r,eps);
% scatter3(r(K,1), r(K,2), r(K,3), 80, 'b', 'filled')
for i=1:length(n)
%     scatter3(r(n(i),1), r(n(i),2), r(n(i),3), 80, 'g', 'filled')
        plot3([r(K,1) r(n(i),1)],...
              [r(K,2) r(n(i),2)],...
              [r(K,3) r(n(i),3)],...
              'LineWidth', 1, 'color', 'black');
end
end
view(3)
%}
end
