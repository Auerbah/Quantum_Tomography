function K = generate_K(lambda, t, optional_precisely)
%������� ���������� ����� ������� K (m, 1), ������������� ������������� ��������,
%������� ����� �� ���������� ��� ���������
%���������:
%   lambda (m, 1) - ��������� �����������
%   t (m, 1) - ����� ���������
%������������ ���������:
%   optional_precisely(=0) - ������������� ��������
%   optional_precisely = 1 -  �� K(k) = lambda(k)*t(k)

if nargin <= 2
	precisely = 0;
else
	precisely = optional_precisely;
end

    if precisely == 0
        K = poissrnd(lambda.*t);
    else
    	K = lambda.*t;
    end
end

