function dm_out = E_a_r_and_p_r_dm(dm_in, T, T1, T2, flag_a_r, flag_p_r)
    s = length(dm_in);
    Ea0 = [[1 0]
          [0 exp(-T/(2*T1))]];
    Ea1 = [[0 sqrt(1-exp(-T/T1))]
          [0 0]];  
    Ep0 = [[1 0]
           [0 exp(-T/(2*T2))]];
    Ep1 = [[0 0]
           [0 sqrt(1-exp(-T/T2))]];
    n = log2(s);
    if flag_a_r == 0 && flag_p_r == 0
        dm_out = dm_in;
        return
    end
    if flag_a_r == 1
        dm_out = zeros(s,s);
        Ea = super_kron(Ea0, Ea1, n);
%         test = eye(s);
%         for i=1:s
%             test = test - Ea(:,:,i)'*Ea(:,:,i);
%         end
%         if test ~= zeros(s,s)
%             error('test ~= zeros(s,s)')
%         end
        for i=1:length(Ea)
            dm_out = dm_out + Ea(:,:,i)*dm_in*Ea(:,:,i)';
        end
        dm_in = dm_out;
    end
    if flag_p_r == 1
        dm_out = zeros(s,s);
        Ep = super_kron(Ep0, Ep1, n);
%         test = eye(s);
%         for i=1:s
%             test = test - Ep(:,:,i)'*Ep(:,:,i);
%         end
%         if test ~= zeros(s,s)
%             error('test ~= zeros(s,s)')
%         end
        for i=1:length(Ep)
            dm_out = dm_out + Ep(:,:,i)*dm_in*Ep(:,:,i)';
        end
    end 
end

