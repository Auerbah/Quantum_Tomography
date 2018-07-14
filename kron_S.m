function N = kron_S(M, S)
%Функция возвращает тензорное произведение матрицы x s раз:
%N =  M x M x M (S = 3)
    N = 1;
    for i=1:S
        N = kron(N,M);
    end
end

