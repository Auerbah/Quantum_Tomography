function Lambda_mix = E_a_r_and_p_r_lambda2(Lambda_in, r, n, T, T1, T2, flag_a_r, flag_p_r)
%Функция возвращает матрицу плотности p_out, полученную в результате
%действия амплитудной релаксации с парметрами t и T1
%Параметры:
%   p_in - изначальная матрица плотности
%   t - время измерения
%   T1 - время амплитудной релаксации
    Ea0 = [[1 0]
          [0 exp(-T/(2*T1))]];
    Ea1 = [[0 sqrt(1-exp(-T/T1))]
          [0 0]];
    Ep0 = [[1 0]
           [0 exp(-T/(2*T2))]];
    Ep1 = [[0 0]
           [0 sqrt(1-exp(-T/T2))]];
    %n - число кубитов, если 2 кубита, то 4 оператора Краусса
    
    if flag_a_r == 0 && flag_p_r == 0
        Lambda_mix = Lambda_in;
        return 
    end
    if flag_a_r == 1
    Ea = super_kron(Ea0, Ea1, n);
    r_s = length(Lambda_in); 
    Ea_mix = zeros(r_s, r_s, length(Ea));
    for i=1:length(Ea)
        Ea_mix(:,:,i) = kron(eye(r),Ea(:,:,i));
    end
    
    Lambda_mix = zeros(length(Lambda_in),length(Lambda_in));
    for i=1:length(Ea)
        Lambda_mix = Lambda_mix + Ea_mix(:,:,i)'*Lambda_in*Ea_mix(:,:,i);
    end
    Lambda_in = Lambda_mix;
    end
    
    if flag_p_r == 1
    Ep = super_kron(Ep0, Ep1, n);
    r_s = length(Lambda_in); 
    Ep_mix = zeros(r_s, r_s, length(Ep));
    for i=1:length(Ep)
        Ep_mix(:,:,i) = kron(eye(r),Ep(:,:,i));
    end
    
    Lambda_mix = zeros(length(Lambda_in),length(Lambda_in));
    for i=1:length(Ep)
        Lambda_mix = Lambda_mix + Ep_mix(:,:,i)'*Lambda_in*Ep_mix(:,:,i);
    end
    end
end

