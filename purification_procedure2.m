function C = purification_procedure2(po)
    
    eps = 1e-7;
    [V, D] = eig(po);
    n = length(po);
    C = [];
    for i = 1:n
        if(abs(D(i,i)) > eps)
            C = [sqrt(D(i,i))*V(:,i); C];
        end
    end
end

