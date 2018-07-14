clc
clear
close
addpath('../Kraus')

%% Параметры шумов
t = 1;
T1 = 1;
T2 = 1;
p_a = 1;

%% Графика
fig = figure('Name','Bloch Sphere','pos',[600 200 600 600]);
figure(fig);
xlabel('$$x$$', 'Interpreter', 'latex');
ylabel('$$y$$', 'Interpreter', 'latex');
zlabel('$$z$$', 'Interpreter', 'latex');
title(sprintf('%s', "Phase relaxation"),'Interpreter', 'none');

%% Начальное состояние
radius = 1;
tet = 3*pi/4;
phi = pi/2;

po = build_po_matrix(radius,tet,phi);
[radius, tet, phi] = return_r_tet_phi_by_po_matrix(po);
r = [radius*cos(phi)*sin(tet),...
     radius*sin(phi)*sin(tet),...
     radius*cos(tet)];

%% Визуализация
type_cmap = 'white';
cmap = colormap(type_cmap);

% axis tight manual % this ensures that getframe() returns a consistent size
filename = 'testAnimated.gif';

t
hold on
[X,Y,Z] = sphere(30);
surf(X, Y, Z, 'EdgeAlpha', 0.1, 'FaceAlpha', 0.1);
scatter3(r(:,1), r(:,2), r(:,3), 80, [0 0 0], 'filled');

az = -90;
el = 0;
view(az, el);

t_time = 0:0.02:T1;
% t_time = t_time.^2;
koeff = 0:1/length(t_time):1;

%% Amplitude and phase relaxation
fprintf("Input density matrix:\n")
po

for i = 1:length(t_time)
      drawnow 
      % Capture the plot as an image 
      frame = getframe(fig); 
      im = frame2im(frame); 
      [imind,cm] = rgb2ind(im,256); 
      % Write to the GIF File 
      if i == 1 
          imwrite(imind,cm,filename,'gif'); 
      else 
          imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0.1); 
      end 
    t = t_time(i);
    k = koeff(i);
    po = E_ampl_and_phase(po, t, T1, p_a, T2);
    [radius, tet, phi] = return_r_tet_phi_by_po_matrix(po);
    r = [radius*cos(phi)*sin(tet),...
         radius*sin(phi)*sin(tet),...
         radius*cos(tet)];
    scatter3(r(:,1), r(:,2), r(:,3), 80, [k 0 0], 'filled');
%     pause(0.1)
end


%% Depolarization noise
% p_time = 0:0.06:0.36;
% koeff = 0:1/length(p_time):1;
% for i = 1:length(p_time)
%     p = p_time(i);
%     k = koeff(i);
%     po = E_depolarization(po, p);
%     [radius, tet, phi] = return_r_tet_phi_by_po_matrix(po);
%     r = [radius*cos(phi)*sin(tet),...
%          radius*sin(phi)*sin(tet),...
%          radius*cos(tet)];
%     scatter3(r(:,1), r(:,2), r(:,3), 80, [k 0 0], 'filled');
%     pause(0.1)
% end

fprintf("Output density matrix:\n")
po

hold off