function pi_val = pi_value(x_exp, x_teor)
n = length(x_exp);
if length(x_teor) ~= n 
    error("lenght(x1) ~= lenght(x2)")
end

chi2 = 0;
for i = 1:n
	chi2 = (x_exp(i) - x_teor(i))^2/x_teor(i);
end
chi2 = chi2*n;
f = @(x) chi2pdf(x, n - 1);
pi_val = integral(f, chi2, inf);
end

