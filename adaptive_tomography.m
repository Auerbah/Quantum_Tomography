% ѕрограмма дл€ восстановлени€ состо€ни€, заданого матрицей плотности
% ¬озможно восстановление как чистого, так и смешанного состо€ни€

tomography_init();
rng(4,'twister');
%% Initial parameters

% Main parameters
    S = 1;                % Number of qubits, if S = 1 (1 qubit)
    x = X_Tetrahedron;    % Chosen hardware matrix
    dm = build_dm(0.9999, pi/4, 5*pi/3); %Chosen quantum state
    n_exp = 1000;          % Number of experiments
    n = 1e6;              % Number of samples in an each experiment
    r = 2;
% Noises
    T = 0;                % Time of measurements
    T1 = 1; flag_a_r = 0; % Time of Ampl relaxation
    T2 = 1; flag_p_r = 0; % Time of Phase relaxation


% Initialize other parameters
[a, counter, dF_set, pi_val_set] = tomography_init2(2^S + 1, n_exp);

%% Main values

% Hardware matrix for chosen dimension
x = kron_S(x, S);

% Number of rows "m" and dimension of system "s" for chosen HW
[m, s] = size(x);
        

% Probabilities of events
lambda = lambda_for_protocol_noise(x, m, dm, T, T1, T2, flag_a_r, flag_p_r);

%% Main cycle


% Show initial and chosen ranks, show m and s
% tomography_info(rank(dm), r, m, s)

alpha0 = [0.001 0.01 0.1 0.25 1/3 0.5 2/3 0.75 0.9 0.99 0. 999];
 
alpha1 = [0.001 0.01 0.1 0.25 1/3 0.5 2/3 0.75 0.9 0.99 0. 999];

DF = zeros(length(alpha0), length(alpha1));
for  i=1:length(alpha0)
for  j=1:length(alpha1)
    
alpha_0 = alpha0(i);
alpha_1 = alpha1(j);
if alpha_0 + alpha_1 > 0
    continue
end
alpha_2 = 1 - (alpha_0 + alpha_1); 
    
% Exposition time for each row of HW
t = time_for_protocol(n * alpha_0, s, m, 'uniform');

for n_e=1:n_exp

% Show progress
progress_info(n_e, 100)

% Zero iteration
    K = generate_K(lambda, t, 0);
    dm_0 = tomography(r, s, x, m, t, K, counter, a, T, T1, T2, flag_a_r, flag_p_r);
    dF_set(n_e, 1) = 1 - fidelity(dm, dm_0);

% First iteration    
    [x_new_1, t_new_1] = new_protocol_by_lorenz_transformation(x, t, dm_0, n*alpha_1);
    lambda_new_1 = lambda_for_protocol_noise(x_new_1, m, dm, T, T1, T2, flag_a_r, flag_p_r);
    K_new_1 = generate_K(lambda_new_1, t_new_1, 0);
    dm_1 = tomography(r, s, [x;x_new_1], 2*m, [t;t_new_1], [K;K_new_1], counter, a, T, T1, T2, flag_a_r, flag_p_r);
    dF_set(n_e, 2) = 1 - fidelity(dm, dm_1);

%Second iteration
    [x_new_2, t_new_2] = new_protocol_by_lorenz_transformation(x, t, dm_1, n*alpha_2);
    lambda_new_2 = lambda_for_protocol_noise(x_new_2, m, dm, T, T1, T2, flag_a_r, flag_p_r);
    K_new_2 = generate_K(lambda_new_2, t_new_2, 0);
    dm_2 = tomography(r, s, [x;x_new_1;x_new_2], 3*m, [t;t_new_1;t_new_2], [K;K_new_1;K_new_2], counter, a, T, T1, T2, flag_a_r, flag_p_r);
    dF_set(n_e, 3) = 1 - fidelity(dm, dm_2);

end

DF(i,j) = mean(dF_set(:,2));
%%
fprintf("\nRESULTS:\n")
fprintf("alpha0   = %0.3f\n", alpha_0)
fprintf("alpha1   = %0.3f\n", alpha_1)
fprintf("alpha2   = %0.3f\n", alpha_2)

fprintf("dF_exp0  = %0.8f\n", mean(dF_set(:,1)))
fprintf("dF_exp1  = %0.8f\n",  mean(dF_set(:,2)))
fprintf("dF_exp2  = %0.8f\n",  mean(dF_set(:,3)))

 end
 end

plot(alpha0, alpha1, DF)

%{
alpha0   = 0.010
alpha1   = 0.100
alpha2   = 0.890
dF_exp0  = 0.00129434
dF_exp1  = 0.00003931
dF_exp2  = 0.00000195

alpha0   = 0.100
alpha1   = 0.100
alpha2   = 0.800
dF_exp0  = 0.00008229
dF_exp1  = 0.00000905
dF_exp2  = 0.00000168

alpha0   = 0.450
alpha1   = 0.100
alpha2   = 0.450
dF_exp0  = 0.00001676
dF_exp1  = 0.00000421
dF_exp2  = 0.00000174

alpha0   = 0.100
alpha1   = 0.450
alpha2   = 0.450
dF_exp0  = 0.00007532
dF_exp1  = 0.00000309
dF_exp2  = 0.00000165

alpha0   = 0.001
alpha1   = 0.010
alpha2   = 0.989
dF_exp0  = 0.00469679
dF_exp1  = 0.00065958
dF_exp2  = 0.00000227

alpha0   = 0.333
alpha1   = 0.333
alpha2   = 0.333
dF_exp0  = 0.00002514
dF_exp1  = 0.00000270
dF_exp2  = 0.00000175
%}