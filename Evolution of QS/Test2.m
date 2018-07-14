clc
clear
close

h = plot(NaN,NaN,'o'); % получаем указатель (handle) точки
axis([0 1000 0 1000]) % задаем границы области графика
for t = 0:0.1:10;
  x=10*t^2;
  y=100*t;
  set(h, 'XData',x, 'YData',y);
  pause(0.01)
end