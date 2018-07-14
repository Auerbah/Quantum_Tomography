function C1 = iteraction_procedure(I, C1, C_tmp, Lambda, K, r, s, m, a, counter)
%Функция возврящает результат итерационной процедуры
%C1 (s_special, 1)
%Параметры:
%   I ()
%   Lambda ()
%   K ()
%   r
%   s
%   m
%   a
%   counter
% p = [];
% Основной цикл
while((1-(abs(C_tmp'*C1))^2)>10^(-10) && counter > 0)
%     p = [p, log((abs(C_tmp'*C1))^2)];
    counter = counter - 1;
    C_tmp = C1;
	L2 = zeros(m,1);

    for k=1:m
        L2(k) = C1' * Lambda(:,:,k) * C1;
    end

	% Матрица, зависящая от С
	J = zeros(r*s, r*s); 
    for k=1:m
        J = J + Lambda(:,:,k) * K(k) / L2(k);
    end
    
	C1 = (1 - a) * I^(-1) * J * C1 + a * C1;
	C1 = C1 / norm(C1);
end
% counter
% plot(p)
end

