function C = purification_procedure3(po)
    
    eps = 1e-7;
    [V, D] = eig(po);
    n = length(po);
    c = [];
    for i = 1:n
        if(abs(D(i,i)) > eps)
            c = [sqrt(D(i,i))*V(:,i) c];
        end
    end
    r = rank(po);
    e = zeros(r,r);
    C = zeros(r*n,1);
    for i=1:r
        e(:,i) = basis(r,i);
    end
    for i=1:r
        C = C + kron(e(:,i),c(:,i));
    end
end