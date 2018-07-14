% ��������� ��� �������������� ���������, �������� �������� ���������
% �������� �������������� ��� �������, ��� � ���������� ���������
clc
clear
close all
addpath('../Protocols');
addpath('Kraus')

rng(1,'twister');

%% ��������� �������
S = 1;             % ����� ������� S = 1 - ���� �����
n_exp = 1000;      % ����� �������������
n = 1000;          % ������� � ������ ������������
a = 0.5;           % �������� ����������
x = X_Octahedron; % ��������� ��������
T = 1;
T1 = 1;

%% ������������� ���������, �� ������ ������� �������� ������� ���������
psi0 = [1; 0];
psi1 = [0; 1];

phi0 = 1 / sqrt(2) * [1; 1];
phi1 = 1 / sqrt(2) * [1; -1];

chi0 = 1 / sqrt(2) * [1; 1i];
chi1 = 1 / sqrt(2) * [1; -1i];

%% ��������� ���������
psi = kron_S(chi0, S)
po = psi*psi'    % ������� ���������

x = kron_S(x, S);

[m, s] = size(x) % m - ����� ����� � ���������� �������
                 % s - ����������� �������
                 
%% ������� ������� ���������� ��� ��������� � ������ ������
t = time_for_protocol(n, s, m, 'uniform');
                 
%% �������������� � ������ r = 1
% ���� ����������� ����� � ��������� �����,
% �� ���� ����� ��������������
r = rank(po);
if(r ~= 1)
    error('rank != 1')
end

%% ������ ����������� ����������� ������� {r = 1}
lambda = lambda_for_protocol_noise(x, m, po, T, T1);

%% ����� ���������� ������� {r = 1}
X = X_for_mixed(r, x);

%% ����� ��������� {r = 1}
Lambda = Lambda_for_mixed_noise(X, T, T1);
sum = zeros(2,2);
for i=1:m
    sum = sum + Lambda(:,:,i);
end
sum

%% ������� I {r = 1}
I = I_for_protocol(Lambda, r, s, t, m);

% ����� ����������� ������ ��������
dF_set = zeros(n_exp,1); 
 
for n_e=1:n_exp
    
if(mod(n_e,100) == 0)
    fprintf('��������� %d �������������\n', n_e);
end
%% ��������� ����� �������, ������������� ������������� ��������
K = generate_K(lambda, t);

%% ������������ ��������� 
% C1 - ������ �����������
C1 = first_approximation(x, t, K, r);
% C1 = psi;

% C_tmp - ��������� ����������
C_tmp = C_temporary(r, s);
%����� �����
counter = 3;
C1 = iteraction_procedure(I, C1, C_tmp, Lambda, K, r, s, m, a, counter);

%% ���������, �������� ��������� ��������
po1 = reverse_purification(C1, r, s);

%% ���������
dF_set(n_e) = 1 - fidelity(po, po1);

end

[dF_average_teor, d] = dF_mixed_with_noise(x, po, t, n, r);
dF_average_exp = mean(dF_set)
dF_average_teor
%% ���������� ��������
plot_hist_and_generalchi2pdf(dF_set, d, n_exp)
