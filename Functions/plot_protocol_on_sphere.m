function plot_protocol_on_sphere(x, size, c)
    for i=1:length(x)
        dm = x(i,:)'*x(i,:);
        plot_dot_on_sphere(dm, size, c)
    end
end

