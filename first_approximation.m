function C = first_approximation( X, t, K, r )

% Нахождеине 1го приближениея для итерационной процедуры
% методом линейной инверсии Мура-Пенроуза

[m, s] = size(X);
B = zeros(m, s*s);
for i=1:m
    B(i,:) = t(i)*kron(conj(X(i,:)),X(i,:));
end

T = K;

[U,S,V] = svd(B);
Q = U'*T;
S = S(1:end-m+s*s,:);
Q = Q(1:end-m+s*s);
f = S^(-1)*Q;
po1 = V*f;
po = [];
for i=1:s
    po = [po, po1((i-1)*s+1:i*s)];
end

% if r == 1
    [U,S,V] = svd(po);
    C = [];
    for i=1:r
        C = [S(i,i)*U(:,i); C];
    end
    C = C/norm(C);
% else
%     C = purification_procedure2(po);
% end

end

