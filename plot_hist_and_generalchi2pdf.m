function plot_hist_and_generalchi2pdf(protocol, dF_set, d, n, n_exp, s, r, T, T1, T2)
%Функция строит гистограмму потерь точности dF = (1-F)
%И теоретическую функцию распределения потерь точности
%Параметры:
%   dF_set (n_exp, 1) - набор экспериментальных данных потерь
%   d (m_special, 1) - набор коэффициентов, полученных исходя из матрицы информации
    dF_set = sort(dF_set);
    min_x = min(dF_set);
    max_x = max(dF_set);
    d_x = (max_x - min_x)/100;
    x_p = min_x:d_x:max_x;
    p = chi2pdf_general_bogdanov(x_p,d);
    fig = figure('Name','Tomography','pos',[700 200 700 600]);
    figure(fig);
    h = histogram(dF_set,25,'FaceColor', [30/255,144/255,255/255]);
    p = p * 4 * n_exp / sum(p);
%     txt1 = '\rho = 0.5\mid0\rangle\langle0\mid+0.5\mid1\rangle\langle1\mid';   
%     tx1 = text(max(dF_set)*0.9, max(p)*0.7, txt1, 'FontSize',14,'Interpreter','tex','HorizontalAlignment','right'); 
%     txt2 = sprintf("\nn = %d, N_{exp} = %d\nt = %0.2f, T_{1} = %0.2f, T_{2} = %0.2f", n, n_exp, T, T1, T2);   
%     tx2 = text(max(dF_set)*0.5, max(p)*0.5, txt2, 'FontSize',14,'Interpreter','tex','HorizontalAlignment','left'); 
%     title(protocol + " protocol", "FontSize", 14)
%     xlabel('$$1-F$$', 'Interpreter', 'latex');
%     ylabel('$$P$$', 'Interpreter', 'latex');
%     fprintf("\nsum(p) = 2*n_exp = %.0f\n",sum(p))
%     hold on
%     plot(x_p, p, 'r', 'LineWidth', 3, 'Color', [178/255, 34/255, 34/255])
%     hold off
%     legend({'Численный эксперимент','Теория'},'FontSize',14,'TextColor','black')
%     sum(h.Values)
    %% Критерий согласия Пирсона
    %{
    n_i = h.Values / n_exp;
    p_i = zeros(1,h.NumBins);
%     h.NumBins
    for i=1:h.NumBins
        x_start = h.BinEdges(i);
        x_end = h.BinEdges(i+1);
        d_x_tmp = (x_end - x_start)/1000;
        x_tmp = x_start:d_x_tmp:x_end;
        y_tmp = chi2pdf_general_bogdanov(x_tmp,d);
        p_i(i) = trapz(x_tmp,y_tmp);
    end
    [res, p, stat] = chi2gof(1:length(p_i), 'Frequency', n_i*n_exp, 'Expected', p_i*n_exp, 'Alpha', 0.05, 'NParams', 1)
%}
    
    createfigure(dF_set, x_p, p)
end


