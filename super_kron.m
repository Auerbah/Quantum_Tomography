function y = super_kron(x1, x2, n)
s = 2^n;
index = 0:1:s-1;
max = strlength(dec2bin(s-1));
arr = strings(1,s);

for i=1:s
    arr(i) = dec2bin(index(i));
  
    while strlength(arr(i)) < max
        arr(i) = '0' + arr(i);
    end
end
arr = char(arr);

y = zeros(s, s, s);
for i=1:s
    buf = arr(:,:,i);
    y_tmp = 1;
    for k=1:max
        if buf(k) == '0'
            y_tmp = kron(y_tmp,x1);
        else
            y_tmp = kron(y_tmp,x2);
        end
    end
    y(:,:,i) = y_tmp;
end


