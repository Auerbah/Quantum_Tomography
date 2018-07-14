clc
clear
close all
addpath('../../Protocols');
addpath('../Kraus');
addpath('../Evolution of QS');
rng(3,'twister');

global d
d = 1e-10;
% ������ ��������� �������
chi0 = 1 / sqrt(2) * [1; 1i];
chi1 = 1 / sqrt(2) * [1; -1i];

p = 0.5;
% ������� ���������
dm = p * chi0 * chi0' + (1 - p) * chi1 * chi1';

x = X_Tetrahedron;
T = 1;             % ����� ���������
T1 = 1;            % ����� ����������� ����������
n = 1000;          % ������� � ������ ������������

[m, s] = size(x); % m - ����� ����� � ���������� �������
                 % s - ����������� �������

t = time_for_protocol(n, s, m, 'uniform');
lambda = lambda_for_protocol_noise(x, m, dm, T, T1);
K = generate_K(lambda, t, 1);

% ������������, ��� �������������� ������ ����� � 
% ������� �����. � ������� LHE ��� ����� ������
dm_1 = chi1 * chi1';
LHE1 = likelyhood_log_mixed(dm_1, x, t, K, T, T1);
[radius, tet, phi] = return_r_tet_phi_by_po_matrix(dm_1);
fprintf("Initial:\n")
fprintf("r = %0.10f, p = %0.5f\nLHE = %0.10f\n",radius,0,LHE1)

% ������ ���������, ������������� ���������, �� ���� ��� �����
% ���������, ���������� ����� � ������� �����. ��� ���
% ����������, � ���. ��� ��� ��� � ��� ��� ��������� �����
[phi, tet] = return_phi_tet_2(chi1);
psi_2 = psi_phi_tet(phi,tet);
equals_psi(chi1, psi_2)
psi_1 = psi_phi_tet(pi+phi,pi-tet);
equals_psi(chi0, psi_1)
if abs(psi_1'*psi_2)^2 > d
    error("��������� �� ������������")
end

% ������� ������ ��� ������� ���������
Lambda = Lambda_for_mixed_noise3(x, 1, T, T1);
% � ����������� ��������� ���� ������:
l_1 = zeros(m,1);
l_2 = zeros(m,1);
for i = 1:m
    l_1(i) = real(psi_1'*Lambda(:,:,i)*psi_1);
    l_2(i) = real(psi_2'*Lambda(:,:,i)*psi_2);
end
% �������, ��� ����� ������ ������� ������� ������� ������ ����
p_new = 0;
fprintf("\n���� �������:\n");

for j=1:10
    p_new = p_new - foo(K, p_new, l_1, l_2, t)/grad_foo2(K, p_new, l_1, l_2, t);
    fprintf("p = %0.10f\n",p_new);
end
fprintf("p = %0.10f\n",p_new);
% � ����� �������� �����-�� dm_2
dm_2 = p_new * psi_1 * psi_1' + (1 - p_new) * psi_2 * psi_2';
fprintf("\nMy p:\n")
LHE = likelyhood_log_mixed(dm_2, x, t, K, T, T1);
[radius, tet, phi] = return_r_tet_phi_by_po_matrix(dm_2);
fprintf("r = %0.10f, p = %0.5f\nLHE = %0.10f\n",radius,p_new,LHE)
  

% � �������� �� � ����������:
fprintf("\nOriginal p:\n")
LHE = likelyhood_log_mixed(dm, x, t, K, T, T1);
[radius, tet, phi] = return_r_tet_phi_by_po_matrix(dm);
fprintf("r = %0.10f, p = %0.5f\nLHE = %0.10f\n",radius,p,LHE)
  

% ��������� �������
dm_3 = additional_procedure(dm_1, x, t, K, T, T1);


function f = foo(K, p, l_1, l_2, t)
    m = length(K);
    f = 0;
    for i = 1:m
        f = f + (K(i)/(p*(l_1(i)-l_2(i))+l_2(i))-t(i))*(l_1(i)-l_2(i));
    end
end

function f = grad_foo(K, p, l_1, l_2, t)
    m = length(K);
    f = 0;
    for i = 1:m
        f = f + K(i)*((l_1(i)-l_2(i))/(p*(l_1(i)-l_2(i))+l_2(i)))^2;
    end
end

function f = grad_foo2(K, p, l_1, l_2, t)
    dp = 0.001;
    f = (foo(K, p+dp, l_1, l_2, t) - foo(K, p, l_1, l_2, t))/dp;
end

% �������������� �������
function equals_psi(psi0, psi1)
    global d
    res = psi0 - psi1;
    for i=1:length(psi0)
    if abs(psi0(i)-psi1(i))>d
        psi0
        psi1
        error('not equal')
    end
    end
end