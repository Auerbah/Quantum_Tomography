% ��������� ��� �������������� ���������, �������� �������� ���������
% �������� �������������� ��� �������, ��� � ���������� ���������
clc
clear
close all
addpath('../Protocols');
addpath('Kraus');
addpath('Evolution of QS');
rng(3,'twister');

%% ��������� �������
n_exp = 1000;      % ����� �������������
n = 10000;         % ������� � ������ ������������
a = 0.5;           % �������� ����������

T = 0.5;             % ����� ���������
T1 = 1;            % ����� ����������� ����������
T2 = 1;            % ����� ������� ����������
flag_a_r = 1;
flag_p_r = 1;
%par : T T1 T2 flag_a_r flag_p_r S r
par = [
       [T T1 T2 flag_a_r flag_p_r 1 1]
%        [T T1 T2 flag_a_r flag_p_r 1 2]
       [T T1 T2 flag_a_r flag_p_r 2 1]
%        [T T1 T2 flag_a_r flag_p_r 2 2]
%        [T T1 T2 flag_a_r flag_p_r 2 3]
%        [T T1 T2 flag_a_r flag_p_r 2 4]
       [T T1 T2 flag_a_r flag_p_r 3 1]
%        [T T1 T2 flag_a_r flag_p_r 3 2]
%        [T T1 T2 flag_a_r flag_p_r 3 4]
%        [T T1 T2 flag_a_r flag_p_r 3 8]
       ];
 n_par = size(par);


%% ��������� ���������



%% ���������� ��������
fig = figure('Name','Tomography','pos',[700 200 700 600]);
figure(fig);
% legend({'Theory'},'FontSize',12,'TextColor','black')
title("Octahedron protocol", "FontSize", 14)
xlabel('$$1-F$$', 'Interpreter', 'latex');
ylabel('$$P$$', 'Interpreter', 'latex');

min_x = min(0);
max_x = max(0.005);
d_x = (max_x - min_x)/100;
x_p = min_x:d_x:max_x;
for i = 1:n_par
    x = X_Octahedron;  % ��������� ��������
    T = par(i,1); T1 = par(i,2); T2 = par(i,3);
    flag_a_r = par(i,4); flag_p_r = par(i,5);
    S = par(i,6); rank_dm = par(i,7);
    dm = create_po_matrix(rank_dm,S); % ������� ���������
    x = kron_S(x, S);
    [m, s] = size(x);
    t = time_for_protocol(n, s, m, 'uniform');
    [dF_average_teor, d] = dF_mixed_with_noise(x, dm, t, T, T1,T2,flag_a_r,flag_p_r, n, rank_dm);
    p = chi2pdf_general_bogdanov(x_p,d);
    p = p * 4 * n_exp / sum(p);
%     txt2 = sprintf("\nn = %d, N_{exp} = %d\nt = %0.2f, T_{1} = %0.2f, T_{2} = %0.2f", n, n_exp, T, T1, T2);   
%     tx2 = text(max_x*0.5, max(p)*0.5, txt2, 'FontSize',14,'Interpreter','tex','HorizontalAlignment','left'); 
%     txt1 = '\leftarrow r = ' + sprintf("%d", rank_dm);
%     y1 = max(p);
%     for j=1:length(p)
%         if p(j) == y1
%             x1 = x_p(j + 2);
%             break
%         end
%     end
    hold on
    if i ~= 1
    plot(x_p, p, 'LineWidth', 3)
    else 
        plot(x_p(2:end), p(2:end), 'LineWidth', 3)
    end
    hold off
%     text(x1,y1,txt1,'FontSize', 14,'Color',[1 0 0])
    fprintf("dF_average_teor = %f\n", dF_average_teor)
end
txt2 = sprintf("\nN = %d, N_{exp} = %d", n, n_exp);   
tx2 = text(2.85*10^-3, 800, txt2, 'FontSize',14,'Interpreter','tex','HorizontalAlignment','left'); 
  
legend({'n = 1, \langle1-F\rangle = 1*10^{-4}',...
        'n = 2, \langle1-F\rangle = 3*10^{-4}',...
        'n = 3, \langle1-F\rangle = 7*10^{-4}'},'FontSize',14,'TextColor','black')
% txt1 = '\rho = 0.5\mid0\rangle\langle0\mid+0.5\mid1\rangle\langle1\mid';   
% tx1 = text(max_x*0.9, max(p)*0.7, txt1, 'FontSize',14,'Interpreter','tex','HorizontalAlignment','right'); 
