function Lambda_noised = Lambda_for_mixed_noise3(x,r,T,T1,T2,flag_a_r,flag_p_r)
% ‘ункуи€ возвравращает набор проекторов дл€ смешанного состо€ни€,
% исход€ из аппаратной матрицы дл€ смешанных состо€ний

% –азмерность такой аппаратной матрицы была:
    [m, s] = size(x);	% m - число строк в исходной аппаратной матрице
                        % s - размерность системы
                        % r - ранг
    X = X_for_mixed(r, x);
    Lambda_pure = zeros(r*s,r*s,m);
    for j=1:m
        for l=1:r
            Lambda_pure(:,:,j) = Lambda_pure(:,:,j)+X(j,:,l)'*X(j,:,l);
        end
    end
    
    Lambda_noised = zeros(r*s,r*s,m);
    for j=1:m
%         Lambda_noised(:,:,j) = E_ampl_lambda2(Lambda_pure(:,:,j),r, log2(s),T,T1);
        Lambda_noised(:,:,j) = E_a_r_and_p_r_lambda2(Lambda_pure(:,:,j),r, log2(s),T,T1,T2,flag_a_r,flag_p_r);
    end
end

