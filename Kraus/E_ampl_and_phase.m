function dm = E_ampl_and_phase(dm, t, T1, T2, flag_a_r, flag_p_r)
%Функция возвращает матрицу плотности p_out, полученную в результате
%действия амплитудной релаксации с парметрами t и T1
%Параметры:
%   p_in - изначальная матрица плотности
%   t - время измерения
%   T1 - время амплитудной релаксации             
    Ea0 = [[1 0]
           [0 exp(-t/(2*T1))]];
    Ea1 = [[0 sqrt(1-exp(-t/T1))]
           [0 0]];
      
    Ep0 = [[1 0]
           [0 exp(-t/(2*T2))]];
    Ep1 = [[0 0]
           [0 sqrt(1-exp(-t/T2))]];
    S = log2(length(dm));
    if S == 1
        if flag_a_r == 1
        dm = Ea0*dm*Ea0' + Ea1*dm*Ea1';
        end
        if flag_p_r == 1
        dm =  Ep0*dm*Ep0' + Ep1*dm*Ep1';
        end
    elseif S == 2  
        error('S>1')
    else
        error('S>2');
    end
end

