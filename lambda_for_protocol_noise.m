function lambda = lambda_for_protocol_noise(x, m, po, T, T1, T2, flag_a_r, flag_p_r)
%���������� ����������� ����������� ������� ��� ���������
%������� ��������� � �������������� ���������
%lambda (m, 1)
%���������:
%   m - ����� ����� � ���������
%   x (m, s) - ��������� ��������
%   po (s, s) - ������� ���������
    [m1, s1] = size(x);
    s = length(po);
    if(m1 ~= m)
        error('�������� m �� ��������� c ������ ����� � X')
    end
    
    if(s1 ~= s)
        error('����������� ������� ��������� �� ��������� � ������������ X')
    end
    
    lambda = zeros(m,1);
    
    %% ������������ ������� ���������
    Lambda = zeros(s,s,m);
    for k=1:m
        Lambda(:,:,k) = x(k,:)'*x(k,:);
    end 
    
    Lambda_T = zeros(s,s,m);
    for k=1:m
        Lambda_T(:,:,k) = E_a_r_and_p_r_lambda(x(k,:)'*x(k,:), T, T1, T2, flag_a_r, flag_p_r);
%         Lambda_T(:,:,k) = E_ampl_lambda(x(k,:)'*x(k,:), T, T1);
   end 
    
    %% ������������ ������� ��������� (����������)
%     sum = zeros(s,s);
%     for k=1:m
%         sum = sum + Lambda_T(:,:,k);
%     end
%     sum;
%     % ������� ���������� ������� (�������� �������)
%     A = abs(sum-m/2*eye(2));
%     assert(all(A(:) <= 1e-10), "������� ������� �� ���������")
%     
%     % ������� ����, ��� ������� ��������� ��������� � ��������
%     for k=1:m
%         A = Lambda_T(:,:,k)^2 - Lambda_T(:,:,k);
%         assert(all(A(:) <= 1e-10), "������� �������� �� ���������")
%     end
    
    %% ������
    for k=1:m
        lambda(k,1) = real(trace(Lambda_T(:,:,k) * po));
    end
    
    %% ������ ������
%     po_out = E_ampl(po,T,T1);
%     lambda2 = zeros(m,1);
%     for k=1:m
%         lambda2(k,1) = real(trace(Lambda(:,:,k) * po_out));
%     end
%     lambda2-lambda
end