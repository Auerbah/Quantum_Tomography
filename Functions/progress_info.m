function progress_info(n_e, d)
%PROGRESSSINFO
    if(mod(n_e, d) == 0)
        fprintf('\t%d experiments completed\n', n_e);
    end
end

