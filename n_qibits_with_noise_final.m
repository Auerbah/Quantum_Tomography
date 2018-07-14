% Программа для восстановления состояния, заданого матрицей плотности
% Возможно восстановление как чистого, так и смешанного состояния

tomography_init();
rng(4,'twister');
%% Initial parameters

% Main parameters
    S = 2;                % Number of qubits, if S = 1 (1 qubit)
    x = X_Octahedron;    % Chosen hardware matrix
%     dm = build_dm(0.99, pi/4, 5*pi/3); %Chosen quantum state
    dm = 1/4*eye(4);
n_exp = 1000;          % Number of experiments
    n = 1e4;              % Number of samples in an each experiment
%     pp = 0.5;
    
% Noises
    T = 0.5;                % Time of measurements
    T1 = 1; flag_a_r = 1; % Time of Ampl relaxation
    T2 = 1; flag_p_r = 1; % Time of Phase relaxation


% Initialize other parameters
[a, counter, dF_set, pi_val_set] = tomography_init2(2^S + 1, n_exp);

%% Main values

% Hardware matrix for chosen dimension
x = kron_S(x, S);

% Number of rows "m" and dimension of system "s" for chosen HW
[m, s] = size(x);
        
% Exposition time for each row of HW
t = time_for_protocol(n, s, m, 'uniform');

% Probabilities of events
lambda = lambda_for_protocol_noise(x, m, dm, T, T1, T2, flag_a_r, flag_p_r);

% test_n_t_lambda(n, t, lambda, 1e-7);
% dm_noised = E_a_r_and_p_r_dm(dm, T, T1, T2, flag_a_r, flag_p_r);

% [x_new, t_new] = new_protocol_by_lorenz_transformation(x, t, dm, n);

% lambda_new = lambda_for_protocol_noise(x_new, m, dm, T, T1, T2, flag_a_r, flag_p_r);

%% Main cycle
for r = 4:4

% Show initial and chosen ranks, show m and s
% tomography_info(rank(dm), r, m, s)

for n_e=1:n_exp
n_e
% Show progress
progress_info(n_e, 100)

% Генерация числа фотонов, подчиняющихся распределению Пуассона
K = generate_K(lambda, t, 0);

% Tomography
dm_0 = tomography(r, s, x, m, t, K, counter, a, T, T1, T2, flag_a_r, flag_p_r);
% dm_0 = tomography(r, s, [x;x], 2*m, [t;t], [K;K], counter, a, T, T1, T2, flag_a_r, flag_p_r);

% Result
dF_set(n_e, r) = 1 - fidelity(dm, dm_0);
% pi_val_set(n_e) = pi_val(K, r, x, m, s, t, dm_0, T, T1, T2, flag_a_r, flag_p_r);

if r == 1

%     dm_new = tomography_for_weak_measurements(dm, dm_0, 100000);

%     L = find_lorenz_transformation(dm);
%     x_new = new_protocol_by_lorenz_transformation(x, L);
%     lambda_new = lambda_for_protocol_noise(x_new, m, dm, T, T1, T2, flag_a_r, flag_p_r);
%     t_new = t;

    K_new = generate_K(lambda_new, t_new, 0);
    dm_2 = tomography(2, s, x_new, m, t_new, K_new, counter, a, T, T1, T2, flag_a_r, flag_p_r);
%     dm_2 = tomography(2, s, [x;x_new], 2*m, [t;t_new], [K;K_new], counter, a, T, T1, T2, flag_a_r, flag_p_r);
 
    dF_set(n_e, s + 1) = 1 - fidelity(dm, dm_2);

end

end

end

[dF_average_teor_1, d_1] = dF_mixed_with_noise_new(x, dm, t,  T, T1,T2,flag_a_r,flag_p_r, n);
% [dF_average_teor_2, d_2] = dF_mixed_with_noise_new(x_new, dm, t_new, T, T1, T2, flag_a_r, flag_p_r, n);
% [dF_average_teor_1, d_1] = dF_mixed_with_noise_new([x;x], dm, [t;t], T, T1, T2, flag_a_r, flag_p_r, 2*n);
% [dF_average_teor_2, d_2] = dF_mixed_with_noise_new([x;x_new], dm, [t;t_new], T, T1, T2, flag_a_r, flag_p_r, 2*n);

fprintf("\nRESULTS:\n")
fprintf("dF_exp=r=1  = %0.8f\n", mean(dF_set(:,1)))
fprintf("dF_exp_r=2  = %0.8f\n",  mean(dF_set(:,2)))
fprintf("dF_teor_old = %0.8f\n", dF_average_teor_1)
fprintf("dF_exp_new  = %0.8f\n", mean(dF_set(:,3)))
% fprintf("dF_teor_new = %0.8f\n", dF_average_teor_2)
%%
% fig = figure('Name','Bloch Sphere','pos',[700 200 700 700]);
% figure(fig);
% 
% hold on
%     
% plot_sphere()
% 
% plot_dot_on_sphere(dm, 200, 'g')
% plot_dot_on_sphere(build_orthogonal_dm(dm), 200, 'b')    
% plot_protocol_on_sphere(x, 200, [0 0 0]);
%     
% plot_protocol_on_sphere(x_new, 50, 'r');

%% Theoretical results and etc
fig = figure('Name','Tomography','pos',[700 200 700 600]);
figure(fig);
hold on
plot_results("Octahedron", dF_set(:,s), d_1, n, n_exp, s, r, T, T1, T2)
plot_results("Octahedron", XXX.dF_set(:,s), XXX.d_1, XXX.n, XXX.n_exp, XXX.s, XXX.r, XXX.T, XXX.T1, XXX.T2)

% plot_results("Tetrahedron lorenz", dF_set(:,s+1), d_2, n, n_exp, s, r, T, T1, T2)

% save_results("X_Octahedron", S, rank(dm), r, T, T1, T2, n_exp, n);

% end

d_1
d_2
