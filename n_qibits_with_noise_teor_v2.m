% Программа для восстановления состояния, заданого матрицей плотности
% Возможно восстановление как чистого, так и смешанного состояния

tomography_init();

%% Initial parameters

% Main parameters
    S = 1;                % Number of qubits, if S = 1 (1 qubit)
    dm = create_dm(45, 0, 0.01);  % Chosen quantum state
%     r = 2;                % Chosen rank of reconstraction
    x = X_Octahedron;     % Chosen hardware matrix
    n_exp = 100;         % Number of experiments
    n = 100000;            % Number of samples in an each experiment

% Noises
    T = 0;                % Time of measurements
    T1 = 1; flag_a_r = 0; % Time of Ampl relaxation
    T2 = 1; flag_p_r = 0; % Time of Phase relaxation

% Initialize other parameters
[a, counter, dF_set, pi_val_set] = tomography_init2(n_exp);


dF_set2 = dF_set;
%% Main values

% Hardware matrix for chosen dimension
x = kron_S(x, S);

% Number of rows "m" and dimension of system "s" for chosen HW
[m, s] = size(x);
                
% Exposition time for each row of HW
t = time_for_protocol(n, s, m, 'uniform');

% Probabilities of events
lambda = lambda_for_protocol_noise(x, m, dm, T, T1, T2, flag_a_r, flag_p_r);

test_n_t_lambda(n, t, lambda, 1e-7);

L = find_lorenz_transformation(dm);
  
x_new = new_protocol_by_lorenz_transformation(x, L);
    
lambda_new = lambda_for_protocol_noise(x_new, m, dm, T, T1, T2, flag_a_r, flag_p_r);
    
t_new = t;

%% Main cycle
for r = 1:1
    
% Show initial and chosen ranks, show m and s
tomography_info(rank(dm), r, m, s)

for n_e=1:n_exp
    
% Show progress
progress_info(n_e, 100)

% Генерация числа фотонов, подчиняющихся распределению Пуассона
K = generate_K(lambda, t, 0);

% Tomography
dm_0 = tomography(r, s, x, m, t, K, counter, a, T, T1, T2, flag_a_r, flag_p_r);

% Result
dF_set(n_e) = 1 - fidelity(dm, dm_0);

pi_val_set(n_e) = pi_val(K, r, x, m, s, t, dm_0, T, T1, T2, flag_a_r, flag_p_r);

if r == 1

    
%     dm_new = tomography_for_weak_measurements(dm, dm_0, 100000);

%     L = find_lorenz_transformation(dm);
%   
%     x_new = new_protocol_by_lorenz_transformation(x, L);
%     
%     lambda_new = lambda_for_protocol_noise(x_new, m, dm, T, T1, T2, flag_a_r, flag_p_r);
%     t_new = t;
    K_new = generate_K(lambda_new, t_new, 0);
    dm_2 = tomography(r, s, x_new, m, t_new, K_new, counter, a, T, T1, T2, flag_a_r, flag_p_r);
    fprintf("Old dF: %f\n", 1 - fidelity(dm, dm_0))
    fprintf("New dF: %f\n", 1 - fidelity(dm, dm_2))
    
    dF_set2(n_e) = 1 - fidelity(dm, dm_2)
    

end

end

    hold on
    
    plot_sphere()
    
    plot_dot_on_sphere(dm_0, 100, 'b')
    
    plot_protocol_on_sphere(x, 50, [0 0 0]);
    
    plot_protocol_on_sphere(x_new, 50, 'r');
%% Theoretical results and etc

% dF_average_teor = dF_mixed_with_noise(x, dm, t, T, T1, T2, flag_a_r, flag_p_r, n, rank(dm));

% print_results(dF_average_teor, dF_set)

% plot_results("Octahedron", dF_set, d, n, n_exp, s, r, T, T1, T2)

% save_results("X_Octahedron", S, rank(dm), r, T, T1, T2, n_exp, n);

end
