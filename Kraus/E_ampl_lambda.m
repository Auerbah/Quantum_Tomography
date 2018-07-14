function Lambda_mix = E_ampl_lambda(Lambda_in, t, T1)
%Функция возвращает матрицу плотности p_out, полученную в результате
%действия амплитудной релаксации с парметрами t и T1
%Параметры:
%   p_in - изначальная матрица плотности
%   t - время измерения
%   T1 - время амплитудной релаксации
    s = length(Lambda_in);
    E0 = [[1 0]
          [0 exp(-t/(2*T1))]];
    E1 = [[0 sqrt(1-exp(-t/T1))]
          [0 0]];
    n = log2(s);
    E = super_kron(E0, E1, n);
    test = eye(s);
    for i=1:s
        test = test - E(:,:,i)'*E(:,:,i);
    end
    if test ~= zeros(s,s)
        error('test ~= zeros(s,s)')
    end
    Lambda_mix = zeros(s,s);
    
    for i=1:length(E)
        Lambda_mix = Lambda_mix + E(:,:,i)'*Lambda_in*E(:,:,i);
    end
%     Lambda_mix-Lambda_mix2
end

