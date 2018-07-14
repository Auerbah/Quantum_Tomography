function lambda = lambda_for_protocol_noise2(x, m, po, T, T1)
%���������� ����������� ����������� ������� ��� ���������
%������� ��������� � �������������� ���������
%lambda (m, 1)
%���������:
%   m - ����� ����� � ���������
%   x (m, s) - ��������� ��������
%   po (s, s) - ������� ���������
    [m1, s] = size(x);
    s2 = length(po);
    if(m1 ~= m)
        error('�������� m �� ��������� c ������ ����� � X')
    end
    
    if(s ~= s2)
        error('����������� ������� ��������� �� ��������� � ������������ X')
    end
    
    lambda = zeros(m,1);
    
    %% ������������� ���������� ������� ����� ����
    phi = zeros(m,1);
    theta = zeros(m,1);
    s2 = zeros(m,1);
    for k=1:m
        phi(k) = phase(x(k,2))-phase(x(k,1));
        phi(k) = -phi(k);
        theta(k) = 2*acos(abs(x(k,1)));
        s2(k) = (phase(x(k,1))+phase(x(k,2)))/2;
    end
    
    %% ������������ ������� ���������
    Lambda_T = zeros(s,s,m);
    % ���������� ������� ��� ���������� C
    for k=1:m
        buffer = [1+cos(theta(k)),sin(theta(k))*exp(1i*phi(k))*exp(-T/(2*T1));
                  sin(theta(k))*exp(-1i*phi(k))*exp(-T/(2*T1)),1+cos(theta(k))*(1-2*exp(-T/(2*T1)))];
        buffer = buffer * 1 / 2;
        Lambda_T(:,:,k) = buffer;
    end
    
        %% ������������ ������� ��������� (����������)
    sum = zeros(s,s);
    for k=1:m
        sum = sum + Lambda_T(:,:,k);
    end
    sum;
%     % ������� ���������� ������� (�������� �������)
%     A = abs(sum-m/2*eye(2));
%     assert(all(A(:) <= 1e-10), '������� ������� �� ���������')
%     
%     % ������� ����, ��� ������� ��������� ��������� � ��������
%     for k=1:m
%         A = Lambda_T(:,:,k)^2 - Lambda_T(:,:,k);
%         assert(all(A(:) <= 1e-10), '������� �������� �� ���������")
%     end
    
    %% ������
    
    for k=1:m
        lambda(k,1) = real(trace(Lambda_T(:,:,k) * po));
    end
    po_out = E_ampl(po, T, T1);
    lambda_2 = zeros(m,1);
    for k=1:m
        lambda_2(k,1) = real(trace(x(k,:)'*x(k,:) * po_out));
    end
    lambda
    lambda_2
end