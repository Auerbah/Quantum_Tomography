function x = new_protocol_by_lorenz(x, dm)

[U, S]= eig(dm);
cTrue = U * sqrt(S);
L = cTrue^(-1)/sqrt(length(dm));
x = x * L;
for i=1:length(x)
    x(i,:) = x(i,:) / norm(x(i,:));
end

end

