clear
addpath('Functions');
scanning_init();
[filename_mat, X_matrix] = protocols_for_scanning();

%% Main parameters
    r1 = 0.99;   % Radius of Bloch Sphere
    n = 1e6;      % Number of samples
    T = 0;        % Time of measurement
    T1 = 1;       % Time of amplitude relaxation
    T2 = 2;       % Time of phase relaxation
    flag_a_r = 0; % Turn on/off amplitude relaxation
    flag_p_r = 0; % Turn on/off phase relaxation
    N = 300;       % Number of dots on Bloch sphere
    
    type_cmap = 'jet';
    
%% Main Loop
[X,Y,Z] = sphere(N);

for k=1:length(X_matrix)
tic
x = cell2mat(X_matrix(k));
dm = build_dm(r1, pi/4, 5*pi/3);
[m, s] = size(x);
t = time_for_protocol(n, s, m, 'uniform');

[x, t] = new_protocol_by_lorenz_transformation(x, t, dm, n);
% [x, t] = new_protocol2(p, x, dm, n);
Losses = zeros(N + 1, N + 1);

for ii = 1:size(X,1)
    for jj = 1:size(X,2)
        dm_x = build_dm_n(r1*[X(ii,jj),Y(ii,jj),Z(ii,jj)]);
%         Losses(ii,jj) =  n * real(dF_mixed_with_noise_new(x, dm_x, t,  T, T1,T2,flag_a_r,flag_p_r, n));
        Losses(ii,jj) =  real(efficiency(x, dm_x, t, T, T1, T2, flag_a_r, flag_p_r, n));
    end
end

%% Plot Bloch Sphere
fig = figure('Name','Bloch Sphere','pos',[800 200 700 700]);
figure(fig);
hold on
%     plot_bloch_sphere(char(filename_mat(k)), Losses, X, Y, Z, r1, type_cmap)
    plot_bloch_sphere('Tetrahedron', Losses, X, Y, Z, r1, type_cmap)
    plot_protocol_on_sphere(cell2mat(X_matrix(k)), 100, 'r');
    plot_protocol_on_sphere(x, 50, 'y');
        if r1 < 1
        surf(X,Y,Z,'EdgeAlpha',0.00, 'FaceAlpha',0.05);
    end
    view(3)
hold off


%% Print results
print_scanning_results(char(filename_mat(k)), Losses)

end
