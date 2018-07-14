function flag = check_po_matrix( po )
%check_po_matrix ������� ���������� true, ���� ��� ������� ���������

flag = true;

% �������� �����
if abs(trace(po) - 1) > 1e-10
    flag = false;
    fprintf('Trace != 1')
end

% �������� �����������
if po ~= po'
    flag = false;
    fprintf('�� ��������')
end

% �������� �� ������������� ��������������
if all(eig(po) >= -1e-10)
else
    flag = false;
    fprintf('�� �����. ������������')
end

end

