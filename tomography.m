function dm = tomography(r, s, x, m, t, K, counter, a, T, T1, T2, flag_a_r, flag_p_r)
%tomography returns reconstracted quantum state input parameters 
%   r - rank of reconstracting state
%   s - dim of quantum state
%   x - hardware matrix (HM)
%   m - number of rows in HM
%   t - time of exposition
%   K - number of clicks of detectors
%   counter - number of steps in Iteration procedure

%% ���������� ���������� ���������
% C1 - ������ �����������
C1 = first_approximation(x, t, K, r);
% C_tmp - ��������� ����������
C_tmp = C_temporary(r, s);
% ���������
Lambda = Lambda_for_mixed_noise3(x, r, T, T1, T2, flag_a_r, flag_p_r);
% ������� I
I = I_for_protocol(Lambda, r, s, t, m);
% ��������� ������ ���������, ���������� � ���������� ������������
% ���������
C1 = iteraction_procedure(I, C1, C_tmp, Lambda, K, r, s, m, a, counter);
% ���������, �������� ��������� ��������
dm = reverse_purification(C1, r, s);

end

