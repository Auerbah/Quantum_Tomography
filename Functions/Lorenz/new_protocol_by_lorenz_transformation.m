function [x_new, t_new] = new_protocol_by_lorenz_transformation(x, t, dm, n)
    
    L = find_lorenz_transformation(dm);
    x_new = x * L;
    t_new = t;
    for i=1:length(x)
        t_new(i) = t_new(i) * norm(x_new(i,:))^2;  
%         t_new(i) = t_new(i) / norm(x_new(i,:))^2;
        x_new(i,:) = x_new(i,:) / norm(x_new(i,:));
    end
    t_new = t_new * 2*n/sum(t_new);
   
end