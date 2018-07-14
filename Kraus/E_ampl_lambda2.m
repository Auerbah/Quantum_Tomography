function Lambda_mix = E_ampl_lambda2(Lambda_in, r, n, T, T1)
%Функция возвращает матрицу плотности p_out, полученную в результате
%действия амплитудной релаксации с парметрами t и T1
%Параметры:
%   p_in - изначальная матрица плотности
%   t - время измерения
%   T1 - время амплитудной релаксации
    E0 = [[1 0]
          [0 exp(-T/(2*T1))]];
    E1 = [[0 sqrt(1-exp(-T/T1))]
          [0 0]];
    %n - число кубитов, если 2 кубита, то 4 оператора Краусса
    E = super_kron(E0, E1, n);

    r_s = length(Lambda_in);
    
    E_mix = zeros(r_s, r_s, length(E));
    for i=1:length(E)
        E_mix(:,:,i) = kron(eye(r),E(:,:,i));
    end
    
    Lambda_mix = zeros(length(Lambda_in),length(Lambda_in));
    for i=1:length(E)
        Lambda_mix = Lambda_mix + E_mix(:,:,i)'*Lambda_in*E_mix(:,:,i);
    end
end

