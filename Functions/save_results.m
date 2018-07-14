function save_results(protocol, S, r_dm, r, T, T1, T2, n_exp, n)

filename = sprintf("Results_S=%d,r_in=%d,r_ch=%d,T=%d,T1=%d,T2=%d,n_exp=%d,n=%d,x=%s",...
                    S, r_dm, r, T, T1, T2, n_exp,n, protocol);
save(filename);

fprintf("\nResults save to file:\n%s", filename);

end

