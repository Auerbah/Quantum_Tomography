function Lambda = Lambda_for_mixed_noise2(X,r,T,T1)
% ‘ункуи€ возвравращает набор проекторов дл€ смешанного состо€ни€,
% исход€ из аппаратной матрицы дл€ смешанных состо€ний

% –азмерность такой аппаратной матрицы была:
% X = zeros(m, r*s, r);
[m, r_s] = size(X);	% m - число строк в исходной аппаратной матрице
                    % r_s - ранг*размерность системы

Lambda_pure = zeros(r_s, r_s, m);
for i=1:m
    Lambda_pure(:,:,i) = X(i,:)'*X(i,:);
end

Lambda = zeros(r_s,r_s,m);
for i=1:m
% 	Lambda(:,:,i) = E_ampl_lambda2(Lambda_pure(:,:,i),s,T,T1);
    Lambda(:,:,i) = E_ampl_lambda(Lambda_pure(:,:,i),T,T1);
end

end

