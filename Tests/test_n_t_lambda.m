function test_n_t_lambda(n, t, lambda, e)
    if abs(sum(t.*lambda) - n) > e 
        error('sum ~= n')
    end
end

