function I = I_for_protocol(Lambda, r, s, t, m)
%������� ���������� ������������ ������� I
%I (r*s, r*s)
%���������:
%   Lambda ()
%   r - ���� ������� ���������
%   s - ����������� �������
%   t (m, 1) - ����� ����������
    I = zeros(r*s, r*s);
    for k=1:m
        I = I + Lambda(:,:,k) * t(k);
    end
end

