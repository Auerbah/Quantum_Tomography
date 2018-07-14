function [dF, d] = dF_mixed_with_noise(x, dm, t,  T, T1,T2,flag_a_r,flag_p_r, n)
%% Описание функции
% Функция возвряет потери точности восстановления для смешанных состояний
% x - аппаратная матрица
% po - матрица плотности
% t - время экспозиции
% n - число измерений
r = rank(dm);

[m, s] = size(x); % m - число строк в аппаратной матрице
                  % s - размерность системы

lambda = lambda_for_protocol_noise(x, m, dm, T, T1, T2, flag_a_r, flag_p_r);

%% Процедура очищения
C = purification_procedure2(dm);

% Величина, необходимая для нахождения матрицы информации
c = [real(C); imag(C)];

%% Новые проекторы
Lambda2 = Lambda_for_mixed_noise3(x,r,T,T1,T2,flag_a_r,flag_p_r);

%% Формирование матрицы информации
% Матрица информации
H = zeros(2*r*s,2*r*s);

Lambda = zeros(2*r*s,2*r*s,m);
for k = 1:m
    Lambda(:,:,k) = [[real(Lambda2(:,:,k)) -imag(Lambda2(:,:,k))]
                     [imag(Lambda2(:,:,k))  real(Lambda2(:,:,k))]];
end

for k=1:m
	if( lambda(k) > 10^(-30))
        H = H + t(k)/lambda(k)*(Lambda(:,:,k)*c)*(Lambda(:,:,k)*c)';
    end
end
H = 2 * H;
% if abs(c'*H*c - 2*n) > 1e-3
% error("c'*H*c ~= 2n")
% end

%% Нахождение исходных параметров
S = eig(H);
S = sort(S);
nu = (2 * s - r) * r;
S = S(length(S)-nu+1:end-1);
if S(1) == 0 
    S = S(2:end);
    fprintf("Еще одно значение убрали\n");
end
d = 1./(2.*S);
d = real(d);

%% Результат
dF = sum(d);
if(isinf(dF))
    dF
end
