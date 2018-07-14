function pi_val = pi_val(K, r, x, m, s, t, dm1, T, T1, T2, flag_a_r, flag_p_r)

K2 = lambda_for_protocol_noise(x, m, dm1, T, T1, T2, flag_a_r, flag_p_r).*t;
[res, pi_val] = chi2gof(1:m, 'Frequency', K, 'Expected', K2, 'Alpha', 0.05, 'NParams', (2*s-r)*r);
fprintf("p-val=%f\n", pi_val)
if res ~= 0 
    fprintf("Bad Hypotesis\n")
end
% sum((K-K2).^2)
end

