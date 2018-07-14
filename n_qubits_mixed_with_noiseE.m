% ��������� ��� �������������� ���������, �������� �������� ���������
% �������� �������������� ��� �������, ��� � ���������� ���������
clc
clear
close all
addpath('../Protocols');
addpath('Kraus');
addpath('Evolution of QS');
rng(4,'twister');

%% ��������� �������
S = 1;             % ����� ������� S = 1 - ���� �����
n_exp = 1000;      % ����� �������������
n = 10000;         % ������� � ������ ������������
a = 0.5;           % �������� ����������
x = X_Octahedron;  % ��������� ��������
% x = X_Tetrahedron;
T = 1;             % ����� ���������
T1 = 1;            % ����� ����������� ����������
T2 = 1;            % ����� ������� ����������
flag_a_r = 1;
flag_p_r = 1;
rank_po = 2;

%% ��������� ���������
% po = create_po_matrix(rank_po,S); % ������� ���������
% po = build_po_matrix(0.5, pi/4, pi/4);
phi0 = 1 / sqrt(2) * [1; 1];
phi1 = 1 / sqrt(2) * [1; -1];

p1 = 0.05;
po = p1*phi0*phi0' + (1-p1)*phi1*phi1';
po = 1/2*eye(2);
x = kron_S(x, S);

[m, s] = size(x); % m - ����� ����� � ���������� �������
                  % s - ����������� �������
fprintf("m = %d, s = %d\n", m, s);
                 
%% ������� ������� ���������� ��� ��������� � ������ ������
t = time_for_protocol(n, s, m, 'uniform');
                 
%% �������������� � ������ r = 1
% ���� ����������� ����� � ��������� �����,
% �� ���� ����� ��������������
r = rank(po);
fprintf("�������� ���� �������:\nr = %d\n", r);
if r ~= rank_po 
    error('r ~= rank_po')
end

%% ������ ����������� ����������� ������� {r = 1}
r = 2;
fprintf("��������� ���� �������:\nr = %d\n", r);
lambda = lambda_for_protocol_noise(x, m, po, T, T1, T2, flag_a_r, flag_p_r);

%% ����� ��������� {r = 1}
Lambda = Lambda_for_mixed_noise3(x,r,T,T1,T2,flag_a_r,flag_p_r);

%% ������� I {r = 1}
I = I_for_protocol(Lambda, r, s, t, m);
sum = 0;
for i=1:m
    sum = sum +t(i)*lambda(i);
end
if abs(sum - n) > 1e-7 
    error('sum ~= n')
end

% ����� ����������� ������ ��������
dF_set = zeros(n_exp,1);
pi_val_set = zeros(n_exp,1); 

% ����� ���������� ������� ����������� �� ���������

dF_first_app = zeros(n_exp,1); 

for n_e=1:n_exp
    
if(mod(n_e,100) == 0)
    fprintf('��������� %d �������������\n', n_e);
end
%% ��������� ����� �������, ������������� ������������� ��������
K = generate_K(lambda, t);

%% ������������ ��������� 
% C1 - ������ �����������
C1 = first_approximation(x, t, K, r);
% C_tmp - ��������� ����������
C_tmp = C_temporary(r, s);
%����� �����
counter = 1000;
C1 = iteraction_procedure(I, C1, C_tmp, Lambda, K, r, s, m, a, counter);

%% ���������, �������� ��������� ��������
po1 = reverse_purification(C1, r, s);

%% �������������� ���������
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

%% ���������
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
% ���������� ��������
plot_hist_and_generalchi2pdf("Octahedron", dF_set, d, n, n_exp, s, r, T, T1, T2)

filename = sprintf("Results/S=%d,r=%d,r=%d,T=%d,T1=%d,T2=%d,n_exp=%d,n=%d,x=X_Octahedron",...
                    S,r,r,T,T1,T2,n_exp,n);
save(filename)
% figure
% hist(pi_val_set, 10)

%%
