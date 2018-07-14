function [a, counter, dF_set, pi_val_set] = tomography_init2(s1, n_exp)
a = 0.5;                      % Parameter of convergence
counter = 1000;               % Number of steps in iteraction procedure  
dF_set = zeros(n_exp,3);      % Set of fidelity losses
pi_val_set = zeros(n_exp,s1);  % Set of pi values
end

