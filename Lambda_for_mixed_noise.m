function Lambda = Lambda_for_mixed_noise(X,T,T1)
% ‘ункуи€ возвравращает набор проекторов дл€ смешанного состо€ни€,
% исход€ из аппаратной матрицы дл€ смешанных состо€ний

% –азмерность такой аппаратной матрицы была:
% X = zeros(m, r*s, r);
    [m, r_s, r] = size(X);	% m - число строк в исходной аппаратной матрице
                            % r_s - ранг*размерность системы
                            % r - ранг

    Lambda = zeros(r_s,r_s,m);
    for j=1:m
        for l=1:r
            Lambda(:,:,j) = Lambda(:,:,j)+X(j,:,l)'*X(j,:,l);
        end
    end
    %1ый вариант:
    Labda_T = zeros(r_s,r_s,m);
    for j=1:m
        Lambda(:,:,j) = E_ampl_lambda(Lambda(:,:,j),T,T1);
    end
    %2ой вариант:
%     Lambda_T = zeros(r_s,r_s,m);
%     for j=1:m
%         for l=1:r
%             Lambda_T(:,:,j) = Lambda_T(:,:,j)+E_ampl(X(j,:,l)'*X(j,:,l),T,T1);
%         end
%     end
end

