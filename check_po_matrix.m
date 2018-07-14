function flag = check_po_matrix( po )
%check_po_matrix функция возвращает true, если это матрица плотности

flag = true;

% Проверка следа
if abs(trace(po) - 1) > 1e-10
    flag = false;
    fprintf('Trace != 1')
end

% Проверка эрмитовости
if po ~= po'
    flag = false;
    fprintf('Не эрмитова')
end

% Проверка на положительную определенность
if all(eig(po) >= -1e-10)
else
    flag = false;
    fprintf('Не полож. определенная')
end

end

