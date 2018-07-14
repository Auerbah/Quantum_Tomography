function [dF, d] = dF_mixed_no_noise_new(x, po, t, n, optional_r)
%% Описание функции
% Функция возвряет потери точности восстановления для смешанных состояний
% x - аппаратная матрица
% po - матрица плотности
% t - время экспозиции
% n - число измерений
% 

%% Ранг матрицы 1 <= r <= s
if nargin > 4
	r = optional_r;
else
	r1 = sqrt(2*(trace(po^2)-1/2));
    r = rank(po); % если r = 1 - чистое состояние, r = 2 - смесь полного ранга
    if((0.95 < r1) && (r1 <= 1))
        r = 1;
    end   
end

[m, s] = size(x); % m - число строк в аппаратной матрице
                 % s - размерность системы

%% Процедура очищения
% C = purification_procedure(po, r);
C = purification_procedure2(po);
%% Новая аппаратная матрица
X = X_for_mixed(r, x);

%% Новые проекторы
Lambda1 = Lambda_for_mixed(X);

%% Вероятности
lambda = zeros(m,1);
for k=1:m
    lambda(k,1) = C' * Lambda1(:,:,k) * C;
end

%% Вычисление некоторых величин
% m - число строк в аппаратной матрице, 
% N0 - размерность пространства очищенного вектора
% s - исходная размерность пространства
% r - ранг
% X = zeros(m, r*s, r);
[m, r_s, r] = size(X);
s = fix(r_s/r);

% Величина, необхадимая для нахождения матрицы информации
c = [real(C); imag(C)];

% %% Представление аппаратной матрицы через углы
% phi = zeros(N,1);
% theta = zeros(N,1);
% s = zeros(N,1);
% for k=1:N
%     phi(k) = phase(X(k,2))-phase(X(k,1));
% 	theta(k) = 2*acos(abs(X(k,1)));
%     s(k) = (phase(X(k,1))+phase(X(k,2)))/2;
% end


%% Формирование матрицы информации
% Матрица информации
H = zeros(2*r_s,2*r_s);
s3=2;
% Аппаратная матрица
X2 = zeros(m, s3, 2*r*s, r);
for j=1:m
	for l=1:r
        X2(j,:,:,l) = [real(X(j,:,l)), -imag(X(j,:,l));
                       imag(X(j,:,l)),  real(X(j,:,l))];
    end
end

% Проекторы
Lambda = zeros(2*r*s,2*r*s,m);

for j=1:m
	for l=1:r
        X_tmp = zeros(s,2*r*s);
        for s1=1:s3
            for s2=1:2*r*s
                X_tmp(s1,s2) = X2(j,s1,s2,l);
            end
        end
        Lambda(:,:,j) = Lambda(:,:,j) + X_tmp'*X_tmp;
    end
end

for k=1:m
   if( lambda(k) > 10^(-30))
         H = H + t(k)/lambda(k)*(Lambda(:,:,k)*c)*(Lambda(:,:,k)*c)';
    end
end

H = 2 * H;

%% Нахождение исходных параметров
S = eig(H);
S = sort(S);
nu = (2*s-r)*r;
S = S(length(S)-nu+1:end-1);
d = 1./(2.*S);
d = real(d);
% d = zeros(nu, 1);
% for k=1:length(d)
%    d(k) = 1 / (2 * S(k));
% end

%% Результат
dF = sum(d);
if(isinf(dF))
    dF
end
