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
x = X_Octahedron;  % Выбранный протокол
% x = X_Tetrahedron;
T = 1;             % Время измерений
T1 = 1;            % Время амплитудной релаксации
T2 = 1;            % Время фазовой релаксации
flag_a_r = 1;
flag_p_r = 1;
rank_po = 2;

%% Выбранное состояние
% po = create_po_matrix(rank_po,S); % Матрица плотности
% po = build_po_matrix(0.5, pi/4, pi/4);
phi0 = 1 / sqrt(2) * [1; 1];
phi1 = 1 / sqrt(2) * [1; -1];

p1 = 0.05;
po = p1*phi0*phi0' + (1-p1)*phi1*phi1';
po = 1/2*eye(2);
x = kron_S(x, S);

[m, s] = size(x); % m - число строк в аппаратной матрице
                  % s - размерность системы
fprintf("m = %d, s = %d\n", m, s);
                 
%% Задание времени экспозиции для измерения в каждом базисе
t = time_for_protocol(n, s, m, 'uniform');
                 
%% Восстановление с рангом r = 1
% Если присутсвует смесь с маленьким весом,
% то ранг нужно переопределить
r = rank(po);
fprintf("Исходный ранг матрицы:\nr = %d\n", r);
if r ~= rank_po 
    error('r ~= rank_po')
end

%% Расчет вероятности регистрации событий {r = 1}
r = 2;
fprintf("Выбранный ранг матрицы:\nr = %d\n", r);
lambda = lambda_for_protocol_noise(x, m, po, T, T1, T2, flag_a_r, flag_p_r);

%% Новые проекторы {r = 1}
Lambda = Lambda_for_mixed_noise3(x,r,T,T1,T2,flag_a_r,flag_p_r);

%% Матрица I {r = 1}
I = I_for_protocol(Lambda, r, s, t, m);
sum = 0;
for i=1:m
    sum = sum +t(i)*lambda(i);
end
if abs(sum - n) > 1e-7 
    error('sum ~= n')
end

% Набор результатов потерь точности
dF_set = zeros(n_exp,1);
pi_val_set = zeros(n_exp,1); 

% Набор отклонений первого приблежения от исходного

dF_first_app = zeros(n_exp,1); 

for n_e=1:n_exp
    
if(mod(n_e,100) == 0)
    fprintf('Проведено %d экспериментов\n', n_e);
end
%% Генерация числа фотонов, подчиняющихся распределению Пуассона
K = generate_K(lambda, t);

%% Итерационная процедура 
% C1 - Первое приближение
C1 = first_approximation(x, t, K, r);
% C_tmp - Временная переменная
C_tmp = C_temporary(r, s);
%Число шагов
counter = 1000;
C1 = iteraction_procedure(I, C1, C_tmp, Lambda, K, r, s, m, a, counter);

%% Процедура, обратная процедуре очищения
po1 = reverse_purification(C1, r, s);

%% Дополнительная процедура
% if r == 1
%     po1 = additional_procedure(po1, x, t, K, T, T1);
% end
% 
%     0.6506
%     0.3262
%     0.3182
%     0.7050

% x = [22 19 19 21 24 16 20 19 17 23]
% y = 20*ones(1,10)

%% Результат
dF_set(n_e) = 1 - fidelity(po, po1);
% K1 = lambda.*t;
% K2 = lambda_for_protocol_noise(x, m, po1, T, T1, T2, flag_a_r, flag_p_r).*t;
% [res, pi_val] = chi2gof(1:m, 'Frequency', K1, 'Expected', K2, 'Alpha', 0.05, 'NParams', (2*s-r)*r);
% if res ~= 0 
%     fprintf("Bad Hypotesis\n")
% end
% pi_val_set(n_e) = pi_val;
end
%% 
[dF_average_teor, d] = dF_mixed_with_noise_new(x, po, t, T, T1,T2,flag_a_r,flag_p_r, n);
dF_average_exp = mean(dF_set);
fprintf("\ndF_average_exp  = %f\n", dF_average_exp)
% xxx=0.0001:0.0001:10;
% trapz(xxx, xxx.*chi2pdf_general_bogdanov(xxx,d))

fprintf("dF_average_teor = %f\n", dF_average_teor)
% Построение графиков
plot_hist_and_generalchi2pdf("Octahedron", dF_set, d, n, n_exp, s, r, T, T1, T2)

filename = sprintf("Results/S=%d,r=%d,r=%d,T=%d,T1=%d,T2=%d,n_exp=%d,n=%d,x=X_Octahedron",...
                    S,r,r,T,T1,T2,n_exp,n);
save(filename)
% figure
% hist(pi_val_set, 10)

%%
