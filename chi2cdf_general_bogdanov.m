function y = chi2cdf_general_bogdanov(x,d)
y = zeros(length(x), 1);
for i=1:length(x)
x1 = 0:0.00001:x(i);
% % f = @(x2) chi2pdf_general_bogdanov(x2,d);
% % y = integral(f, 0, x);
y(i) = trapz(x1, chi2pdf_general_bogdanov(x1,d));
% y = y(end)
% x
% y = sin(x)
end
end

