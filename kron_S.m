function N = kron_S(M, S)
%������� ���������� ��������� ������������ ������� x s ���:
%N =  M x M x M (S = 3)
    N = 1;
    for i=1:S
        N = kron(N,M);
    end
end

