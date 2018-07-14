function print_results(dF_average_teor, dF_set)
    dF_average_exp = mean(dF_set);
    fprintf("\nResults:\n", dF_average_exp)
    fprintf("\tdF_average_exp  = %f\n", dF_average_exp)
%     xxx=0.0001:0.0001:10;
%     trapz(xxx, xxx.*chi2pdf_general_bogdanov(xxx,d))
    fprintf("\tdF_average_teor = %f\n", dF_average_teor)
end

