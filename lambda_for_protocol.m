function lambda = lambda_for_protocol(x, m, po)
%���������� ����������� ����������� ������� ��� ���������
%������� ��������� � �������������� ���������
%lambda (m, 1)
%���������:
%   m - ����� ����� � ���������
%   x (m, s) - ��������� ��������
%   po (s, s) - ������� ���������
    [m1, s1] = size(x);
    s = size(po);
    if(m1 ~= m)
        error('�������� m �� ��������� c ������ ����� � X')
    end
    
    if(s1 ~= s)
        error('����������� ������� ��������� �� ��������� � ������������ X')
    end
    
    lambda = zeros(m,1);
    for k=1:m
        lambda(k,1) = real(x(k,:) * po * x(k,:)');
    end
end

