function dm = tomography(r, s, x, m, t, K, counter, a, T, T1, T2, flag_a_r, flag_p_r)
%tomography returns reconstracted quantum state input parameters 
%   r - rank of reconstracting state
%   s - dim of quantum state
%   x - hardware matrix (HM)
%   m - number of rows in HM
%   t - time of exposition
%   K - number of clicks of detectors
%   counter - number of steps in Iteration procedure

%% Томография квантового состояния
% C1 - Первое приближение
C1 = first_approximation(x, t, K, r);
% C_tmp - Временная переменная
C_tmp = C_temporary(r, s);
% Проекторы
Lambda = Lambda_for_mixed_noise3(x, r, T, T1, T2, flag_a_r, flag_p_r);
% Матрица I
I = I_for_protocol(Lambda, r, s, t, m);
% Очищенный вектор состояния, полученный в результате итерационной
% процедуры
C1 = iteraction_procedure(I, C1, C_tmp, Lambda, K, r, s, m, a, counter);
% Процедура, обратная процедуре очищения
dm = reverse_purification(C1, r, s);

end

