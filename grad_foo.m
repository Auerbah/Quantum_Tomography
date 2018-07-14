function f = grad_foo(K, p, l_1, l_2, t)
    dp = 0.001;
    f = (foo(K, p+dp, l_1, l_2, t) - foo(K, p, l_1, l_2, t))/dp;
end