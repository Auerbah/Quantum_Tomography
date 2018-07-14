clc
clear
close

h = plot(NaN,NaN,'o'); % �������� ��������� (handle) �����
axis([0 1000 0 1000]) % ������ ������� ������� �������
for t = 0:0.1:10;
  x=10*t^2;
  y=100*t;
  set(h, 'XData',x, 'YData',y);
  pause(0.01)
end