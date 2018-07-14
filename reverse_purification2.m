function po_a = reverse_purification2(C1, r, s)
%Возвращает матрицу плотности на основе очищенного вектора состояния
%po (s, s) - полученная матрица плотности
%   C1 (s_special, 1) - очищенный вектор состояния
%   r - ранг матрицы плотности
%   s - размерность системы

    po_a = zeros(s,s);
    po = C1*C1';
    for i=1:s
        po_a = po_a + kron(basis(r,i)',eye(s))*po*kron(basis(r,i),eye(s));
    end
%     po_b = zeros(s,s);
%     for i=1:s
%         po_b = po_b + kron(eye(s),basis(r,i)')*po*kron(eye(s),basis(r,i));
%     end
%    po_b
end
