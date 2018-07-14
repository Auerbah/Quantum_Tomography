% Программа для восстановления состояния, заданого матрицей плотности
% Возможно восстановление как чистого, так и смешанного состояния
clc
clear
close all
addpath('../Protocols');
addpath('Kraus');
addpath('Evolution of QS');
rng(4,'twister');

%% Начальные условия
S = 1;             % Число кубитов S = 1 - один кубит
n_exp = 1000;      % Число экспериментов
n = 10000;         % Выборка в каждом эксперименте
a = 0.5;           % Параметр сходимости

x = X_Octahedron;

T = 1;             % Время измерений
T1 = 1;            % Время амплитудной релаксации
T2 = 1;            % Время фазовой релаксации
flag_a_r = 1;
flag_p_r = 1;

dm = 1/2*eye(2);

[m, s] = size(x);
t = time_for_protocol(n, s, m, 'uniform');

%% 

[dF_average_teor_1, d_1] = dF_mixed_with_noise_new(x, dm, t, T, T1, T2, flag_a_r, flag_p_r, n);
fprintf("dF_average_teor_1 = %f\n", dF_average_teor_1)

[dF_average_teor_2, d_2] = dF_mixed_with_noise_new(x, dm, t, T, T1, T2, flag_a_r, flag_p_r, n);
fprintf("dF_average_teor_2 = %f\n", dF_average_teor_2)

[dF_average_teor_3, d_3] = dF_mixed_with_noise_new(x, dm, t, T, T1, T2, flag_a_r, flag_p_r, n);
fprintf("dF_average_teor_3 = %f\n", dF_average_teor_3)

[dF_average_teor_4, d_4] = dF_mixed_with_noise_new(x, dm, t, T, T1, T2, flag_a_r, flag_p_r, n);
fprintf("dF_average_teor_4 = %f\n", dF_average_teor_4)

x_p = 0.1e-3:0.1e-4:10e-3;
p_1 = chi2pdf_general_bogdanov(x_p, d_1);
p_2 = chi2pdf_general_bogdanov(x_p, d_2);
p_3 = chi2pdf_general_bogdanov(x_p, d_3);
p_4 = chi2pdf_general_bogdanov(x_p, d_4);

figure1 = figure('Name','Tomography');

% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');


% Create text
text('Parent',axes1,'HorizontalAlignment','right','FontSize',16,...
    'FontName','Times New Roman',...
    'String','\rho = 0.5\mid0\rangle\langle0\mid+0.5\mid1\rangle\langle1\mid',...
    'Position',[0.00600906136620883 99.0947224413716 0]);

% Create text
text('Parent',axes1,'FontSize',16,'FontName','Times New Roman',...
    'String',{'','n = 10000, N_{exp} = 1000','t = 1.00, T_{1} = 1.00, T_{2} = 1.00'},...
    'Position',[0.00383073541247025 76.9834949885766 0]);

% Create plot
plot(x_p,p_1,'DisplayName','r = 1','HitTest','off','LineWidth',10,...
    'Color',[0.92156862745098 0.396078431372549 0.317647058823529 0.8]);

plot(x_p,p_2,'DisplayName','Протокол куба 6 строк','HitTest','off','LineWidth',10,...
    'Color',[0.92156862745098 0.396078431372549 0.317647058823529 0.8]);

plot(x_p,p_3,'DisplayName','Протокол октаэдра 8 строк','HitTest','off','LineWidth',10,...
    'Color',[0.92156862745098 0.396078431372549 0.317647058823529 0.8]);

% Create ylabel
ylabel('$$P$$','FontName','Times New Roman','Interpreter','latex');

% Create xlabel
xlabel('$$1-F$$','FontName','Times New Roman','Interpreter','latex');

% Create title
title('Сравнение различных протоколов','FontSize',20,'FontName','Times New Roman');

box(axes1,'on');
% Set the remaining axes properties
set(axes1,'Color',[0.913725490196078 0.862745098039216 0.788235294117647],...
    'FontSize',14);
% Create legend
legend1 = legend(axes1,'show');
set(legend1,...
    'Position',[0.717850287907869 0.740972617851471 0.174504153379996 0.138061649871595],...
    'FontSize',14,...
    'EdgeColor',[0.666666666666667 0.666666666666667 0.658823529411765],...
    'Color',[0.666666666666667 0.666666666666667 0.658823529411765]);
