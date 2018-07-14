function f = foo(K, p, l_1, l_2, t)
    m = length(K);
    f = 0;
    for i = 1:m
        f = f + (K(i)/(p*(l_1(i)-l_2(i))+l_2(i))-t(i))*(l_1(i)-l_2(i));
    end
end
