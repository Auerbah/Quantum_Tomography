clc
clear
close
addpath('../Kraus')
addpath('..')

%% Параметры шумов

t_end = 0.5;    %Конечное время измерений
T1 = 1;         %Время амплитуудной релаксации
T2 = 1;         %Время фазовой релаксации

%Наличие шумов: 1 - есть, 0 - нет
flag_a_r = 1;   %Наличие амплитудной релаксации
flag_p_r = 0;   %Наличие фазовой релаксацции

%% Графика
fig = figure('Name','Bloch Sphere','pos',[700 200 800 800]);
figure(fig);
az = 180*3/4;
el = 180/4;
view(az, el);
view(3)
xlabel('$$x$$', 'Interpreter', 'latex');
ylabel('$$y$$', 'Interpreter', 'latex');
zlabel('$$z$$', 'Interpreter', 'latex');
graph_title = "";
if flag_a_r == 1 && flag_p_r == 1
    graph_title = sprintf("Amplitude T_{1} = %d and Phase T_{2} = %d relaxation", T1, T2);
elseif flag_a_r == 1 && flag_p_r == 0
    graph_title = sprintf("Amplitude (T_{1} = %d) relaxation", T1);
elseif flag_a_r == 0 && flag_p_r == 1
    graph_title = sprintf("Phase (T_{2} = %d) relaxation", T2);
end

tit = title(graph_title +  sprintf(" t = %0.2f", 0.00), 'Interpreter', 'tex');
xlim([-1 1])
ylim([-1 1])
zlim([-1 1])

%% Начальное состояние

tet0 = [pi, 3*pi/4, pi/4, 3*pi/4, pi/4, pi/2, pi/2];
phi0 = [0, pi/2, pi/2, 3*pi/2, 3*pi/2, pi/2, 3*pi/2];
tet0 = [tet0, pi, 3*pi/4, pi/4, 3*pi/4, pi/4, pi/2, pi/2];
phi0 = [phi0, 0, 0, 0, pi, pi, 0, pi];
N = length(tet0);
radius0 = ones(1, N);

dm = zeros(2,2,N);
for i=1:N
    dm(:,:,i) = build_po_matrix(radius0(i),tet0(i),phi0(i));
end

radius = zeros(1,N);
tet = zeros(1,N);
phi = zeros(1,N);
r = zeros(N,3);
for i=1:N
    [radius(i), tet(i), phi(i)] = return_r_tet_phi_by_po_matrix(dm(:,:,i));
    
    r(i,:) = [radius(i)*cos(phi(i))*sin(tet(i)),...
              radius(i)*sin(phi(i))*sin(tet(i)),...
              radius(i)*cos(tet(i))];
end



%% Визуализация
type_cmap = 'white';
cmap = colormap(type_cmap);

% axis tight manual % this ensures that getframe() returns a consistent size
filename = 'testAnimated2.gif';

hold on
[X,Y,Z] = sphere(30);
surf(X, Y, Z, 'EdgeAlpha', 0.1, 'FaceAlpha', 0.1);
% txt = texlabel(sprintf("t = %0.2f", 0.00));   
% text_label = text(0.8,0.8,3,txt,'FontSize',14);

for i=1:N
    scatter3(r(i,1), r(i,2), r(i,3), 80, [0 0 0], 'filled')
end
% t_time = [0:0.0001:0.001 0.002:0.001:0.01 0.02:0.01:0.1 0.2:0.1:1]
t_time = 0:0.02:t_end;
% t_time = t_time.^2;
koeff = 0:1/length(t_time):1;

%% Amplitude and phase relaxation

for j = 1:length(t_time)
      az = az + 1;
      view(az, el);
      drawnow 
      % Capture the plot as an image 
      frame = getframe(fig); 
      im = frame2im(frame); 
      [imind,cm] = rgb2ind(im,256); 
      % Write to the GIF File 
      if j == 1 
          imwrite(imind,cm,filename,'gif','Loopcount',inf); 
      else 
          imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0.5); 
      end 
    t = t_time(j);
    k = koeff(j);
    for i=1:N
    dm(:,:,i) = E_a_r_and_p_r_dm(dm(:,:,i), t, T1, T2, flag_a_r, flag_p_r);
%     dm(:,:,i) = rotation_operator( pi/16, 0, 0 )*dm(:,:,i)*...
%                 rotation_operator( pi/16, 0, 0 )';
    
    [radius(i), tet(i), phi(i)] = return_r_tet_phi_by_po_matrix(dm(:,:,i));
 
    r(i,:) = [radius(i)*cos(phi(i))*sin(tet(i)),...
              radius(i)*sin(phi(i))*sin(tet(i)),...
              radius(i)*cos(tet(i))];
    
    tit.String = graph_title +  sprintf(" t = %0.2f", t);
    scatter3(r(i,1), r(i,2), r(i,3), 80, [0 k 0], 'filled');
    end
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

% fprintf("Output density matrix:\n")
% dm

hold off