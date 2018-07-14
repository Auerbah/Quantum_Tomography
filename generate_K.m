function K = generate_K(lambda, t, optional_precisely)
%Функция генерирует число фотонов K (m, 1), подчиняющихся распределению Пуассона,
%которое могло бы получиться при измерении
%Параметры:
%   lambda (m, 1) - амплитуда вероятности
%   t (m, 1) - время измерений
%Опциональные параметры:
%   optional_precisely(=0) - распределение Пуассона
%   optional_precisely = 1 -  то K(k) = lambda(k)*t(k)

if nargin <= 2
	precisely = 0;
else
	precisely = optional_precisely;
end

    if precisely == 0
        K = poissrnd(lambda.*t);
    else
    	K = lambda.*t;
    end
end

