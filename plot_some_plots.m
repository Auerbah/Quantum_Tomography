% Программа для восстановления состояния, заданого матрицей плотности
% Возможно восстановление как чистого, так и смешанного состояния
clc
clear
close all
addpath('../Protocols');
addpath('Kraus');
addpath('Evolution of QS');
rng(3,'twister');

%% Начальные условия
S = 1;             % Число кубитов S = 1 - один кубит
n_exp = 1000;      % Число экспериментов
n = 10000;         % Выборка в каждом эксперименте
a = 0.5;           % Параметр сходимости
x = X_Octahedron;  % Выбранный протокол
% x = X_Tetrahedron;
T = 1;             % Время измерений
T1 = 1;            % Время амплитудной релаксации
T2 = 1;            % Время фазовой релаксации
flag_a_r = 1;
flag_p_r = 1;
rank_po = 2;
x = kron_S(x, S);
[m, s] = size(x); % m - число строк в аппаратной матрице
                  % s - размерность системы
t = time_for_protocol(n, s, m, 'uniform');

%par : T T1 T2 flag_a_r flag_p_r
par = [[0 1 1 1 1]
       [0.25 1 1 0 1]
       [0.25 1 1 1 0]
       [0.25 1 1 1 1]
       [0.5 1 1 1 1]];
 n_par = size(par);


%% Выбранное состояние
po = create_po_matrix(rank_po,S); % Матрица плотности


%% Построение графиков
fig = figure('Name','Tomography','pos',[700 200 700 600]);
figure(fig);
% legend({'Theory'},'FontSize',12,'TextColor','black')
title("Octahedron protocol", "FontSize", 14)
xlabel('$$1-F$$', 'Interpreter', 'latex');
ylabel('$$P$$', 'Interpreter', 'latex');

min_x = min(0);
max_x = max(0.001);
d_x = (max_x - min_x)/100;
x_p = min_x:d_x:max_x;

for i = 1:n_par
    T = par(i,1); T1 = par(i,2); T2 = par(i,3);
    flag_a_r = par(i,4); flag_p_r = par(i,5);
    [dF_average_teor, d] = dF_mixed_with_noise(x, po, t, T, T1,T2,flag_a_r,flag_p_r, n, rank(po));
    p = chi2pdf_general_bogdanov(x_p,d);
    p = p * 4 * n_exp / sum(p);
%     txt2 = sprintf("\nn = %d, N_{exp} = %d\nt = %0.2f, T_{1} = %0.2f, T_{2} = %0.2f", n, n_exp, T, T1, T2);   
%     tx2 = text(max_x*0.5, max(p)*0.5, txt2, 'FontSize',14,'Interpreter','tex','HorizontalAlignment','left'); 
    hold on
    plot(x_p, p, 'LineWidth', 3)
    hold off
    fprintf("dF_average_teor = %f\n", dF_average_teor)
end

txt1 = '\rho = 0.5\mid0\rangle\langle0\mid+0.5\mid1\rangle\langle1\mid';   
tx1 = text(max_x*0.9, max(p)*0.7, txt1, 'FontSize',14,'Interpreter','tex','HorizontalAlignment','right'); 
