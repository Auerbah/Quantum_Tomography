function I = I_for_protocol(Lambda, r, s, t, m)
%функция возвращает неизменяемую матрицу I
%I (r*s, r*s)
%Параметры:
%   Lambda ()
%   r - ранг матрицы плотности
%   s - размерность системы
%   t (m, 1) - время экспозиции
    I = zeros(r*s, r*s);
    for k=1:m
        I = I + Lambda(:,:,k) * t(k);
    end
end

