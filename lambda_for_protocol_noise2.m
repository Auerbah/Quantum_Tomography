function lambda = lambda_for_protocol_noise2(x, m, po, T, T1)
%¬озвращает веро€тности регистрации событий при измерении
%матрицы плотности с использованием протокола
%lambda (m, 1)
%ѕараметры:
%   m - число строк в протоколе
%   x (m, s) - выбранный протокол
%   po (s, s) - матрица плотности
    [m1, s] = size(x);
    s2 = length(po);
    if(m1 ~= m)
        error('Ќеверное m не совпадает c числом строк в X')
    end
    
    if(s ~= s2)
        error('–азмерность матрицы плотности не совпадает с размерностью X')
    end
    
    lambda = zeros(m,1);
    
    %% ѕредставление аппаратной матрицы через углы
    phi = zeros(m,1);
    theta = zeros(m,1);
    s2 = zeros(m,1);
    for k=1:m
        phi(k) = phase(x(k,2))-phase(x(k,1));
        phi(k) = -phi(k);
        theta(k) = 2*acos(abs(x(k,1)));
        s2(k) = (phase(x(k,1))+phase(x(k,2)))/2;
    end
    
    %% ‘ормирование матрицы измерений
    Lambda_T = zeros(s,s,m);
    % Ќеизменна€ матрица дл€ нахождени€ C
    for k=1:m
        buffer = [1+cos(theta(k)),sin(theta(k))*exp(1i*phi(k))*exp(-T/(2*T1));
                  sin(theta(k))*exp(-1i*phi(k))*exp(-T/(2*T1)),1+cos(theta(k))*(1-2*exp(-T/(2*T1)))];
        buffer = buffer * 1 / 2;
        Lambda_T(:,:,k) = buffer;
    end
    
        %% “естирование матрицы измерений (проекторов)
    sum = zeros(s,s);
    for k=1:m
        sum = sum + Lambda_T(:,:,k);
    end
    sum;
%     % ”словие разложение единицы (условние полноты)
%     A = abs(sum-m/2*eye(2));
%     assert(all(A(:) <= 1e-10), 'условие полноты не выполнено')
%     
%     % ”словие того, что квадрат проектора переходит в проектор
%     for k=1:m
%         A = Lambda_T(:,:,k)^2 - Lambda_T(:,:,k);
%         assert(all(A(:) <= 1e-10), 'условие квадрата не выполнено")
%     end
    
    %% Ћ€мбды
    
    for k=1:m
        lambda(k,1) = real(trace(Lambda_T(:,:,k) * po));
    end
    po_out = E_ampl(po, T, T1);
    lambda_2 = zeros(m,1);
    for k=1:m
        lambda_2(k,1) = real(trace(x(k,:)'*x(k,:) * po_out));
    end
    lambda
    lambda_2
end