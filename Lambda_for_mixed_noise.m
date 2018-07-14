function Lambda = Lambda_for_mixed_noise(X,T,T1)
% ������� ������������� ����� ���������� ��� ���������� ���������,
% ������ �� ���������� ������� ��� ��������� ���������

% ����������� ����� ���������� ������� ����:
% X = zeros(m, r*s, r);
    [m, r_s, r] = size(X);	% m - ����� ����� � �������� ���������� �������
                            % r_s - ����*����������� �������
                            % r - ����

    Lambda = zeros(r_s,r_s,m);
    for j=1:m
        for l=1:r
            Lambda(:,:,j) = Lambda(:,:,j)+X(j,:,l)'*X(j,:,l);
        end
    end
    %1�� �������:
    Labda_T = zeros(r_s,r_s,m);
    for j=1:m
        Lambda(:,:,j) = E_ampl_lambda(Lambda(:,:,j),T,T1);
    end
    %2�� �������:
%     Lambda_T = zeros(r_s,r_s,m);
%     for j=1:m
%         for l=1:r
%             Lambda_T(:,:,j) = Lambda_T(:,:,j)+E_ampl(X(j,:,l)'*X(j,:,l),T,T1);
%         end
%     end
end

