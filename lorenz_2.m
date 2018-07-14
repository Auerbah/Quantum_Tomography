clc
clear
close
addpath('Functions/Density Matrix')
addpath('Functions/Lorenz')
addpath('Kraus')
addpath('Functions')
addpath('../Protocols')
x = X_Tetrahedron;
[m, s] = size(x)
n = 1e6;
dm = build_dm(0.99999, pi/4, 5*pi/3);
% dm = [1;0;0;0] + [0;0;0;1] + [0;2i;3i;0]; 
% dm = dm/norm(dm);
% dm = dm*dm';
T=0, T1=1, T2=1, flag_a_r=0, flag_p_r=0

t = time_for_protocol(n, s, m, 'uniform')

[U, S]= eig(dm);
cTrue = U * sqrt(S);
L = cTrue^(-1) / sqrt(2);
L
L*dm*L'

x1 = x * L;
t1 = t;
for i=1:length(x)
    t1(i) = t1(i) * norm(x1(i,:))^2;  
% 	t_new(i) = t_new(i) / norm(x_new(i,:))^2;
	x1(i,:) = x1(i,:) / norm(x1(i,:));
end
t1 = t1 * 2*n/sum(t1);

L = find_lorenz_transformation(dm);
L*dm*L';
L = L/sqrt(trace(L*dm*L'))
L*dm*L'

x2 = x * L;
t2 = t;
for i=1:length(x)
    t2(i) = t2(i) * norm(x2(i,:))^2;  
% 	t_new(i) = t_new(i) / norm(x_new(i,:))^2;
	x2(i,:) = x2(i,:) / norm(x2(i,:));
end
t2 = t2 * 2*n/sum(t2);

fig = figure('Name','Bloch Sphere','pos',[700 200 700 700]);
figure(fig);

hold on
    
plot_sphere()

plot_dot_on_sphere(dm, 200, 'g')
plot_dot_on_sphere(build_orthogonal_dm(dm), 200, 'b')    
plot_protocol_on_sphere(x, 200, [0 0 0]);
    
plot_protocol_on_sphere(x1, 50, 'r');
plot_protocol_on_sphere(x2, 50, 'g');

[dF_average_teor0, d0] = dF_mixed_with_noise_new(x, dm, t, T, T1, T2, flag_a_r, flag_p_r, n);
[dF_average_teor1, d1] = dF_mixed_with_noise_new(x1, dm, t1, T, T1, T2, flag_a_r, flag_p_r, n);
[dF_average_teor2, d2] = dF_mixed_with_noise_new(x2, dm, t2, T, T1, T2, flag_a_r, flag_p_r, n);
dF_average_teor0
dF_average_teor1
dF_average_teor2
