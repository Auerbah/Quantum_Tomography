function dm_2 = additional_procedure(dm, x, t, K, T, T1)

% Смотрим значение логарифма правдоподобия вначале:
LHE_1 = likelyhood_log_mixed(dm, x, t, K, T, T1);

[radius, tet, phi] = return_r_tet_phi_by_po_matrix(dm);
% Строим соcтояние, ортогональное наденному, по идее эту будет
% состояние, отвечающее смеси с меньшим весом. Так что
% убеждаемся, в том. что это так и что эти состояния равны
psi_2 = psi_phi_tet(phi,tet);
psi_1 = psi_phi_tet(pi+phi,pi-tet);
if abs(psi_1'*psi_2)^2 > 1e-10
    error("Состояния не ортогональны")
end

% Находим лямбды для чистого состояния
Lambda = Lambda_for_mixed_noise3(x, 1, T, T1);
% И вероятности измерений этих смесей:
[m, s] = size(x);
l_1 = zeros(m,1);
l_2 = zeros(m,1);
for i = 1:m
    l_1(i) = real(psi_1'*Lambda(:,:,i)*psi_1);
    l_2(i) = real(psi_2'*Lambda(:,:,i)*psi_2);
end
% Говорим, что будем искать решение методом Ньютона вблизи нудя
p_new = 0.5;
for j=1:10
    p_new = p_new - foo(K, p_new, l_1, l_2, t)/grad_foo(K, p_new, l_1, l_2, t);
%     fprintf("p = %0.10f\n",p_new)
end

if p_new < 0 || p_new > 1
%     fprintf('p = %f\n', p_new)
    dm_2 = dm;
    return
end
% В итоге получили какую-то dm_2
dm_2 = p_new * psi_1 * psi_1' + (1 - p_new) * psi_2 * psi_2';
check_po_matrix(dm_2);

LHE_2 = likelyhood_log_mixed(dm_2, x, t, K, T, T1);

if LHE_1 > LHE_2
    error("правдоподобие упало")
end


end

