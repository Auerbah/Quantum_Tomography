function y = likelyhood_log_mixed(po, x, t, K, T, T1)
    y = 0;
    [m, s] = size(x);
    lambda = lambda_for_protocol_noise(x, m, po, T, T1);
    for i = 1:m
        y=y+K(i)*log(lambda(i)*t(i))-lambda(i)*t(i)-(K(i)*log(K(i))-K(i));
    end
end
