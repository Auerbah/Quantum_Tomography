function [dF, d] = dF_mixed_no_noise_new(x, po, t, n, optional_r)
%% �������� �������
% ������� �������� ������ �������� �������������� ��� ��������� ���������
% x - ���������� �������
% po - ������� ���������
% t - ����� ����������
% n - ����� ���������
% 

%% ���� ������� 1 <= r <= s
if nargin > 4
	r = optional_r;
else
	r1 = sqrt(2*(trace(po^2)-1/2));
    r = rank(po); % ���� r = 1 - ������ ���������, r = 2 - ����� ������� �����
    if((0.95 < r1) && (r1 <= 1))
        r = 1;
    end   
end

[m, s] = size(x); % m - ����� ����� � ���������� �������
                 % s - ����������� �������

%% ��������� ��������
% C = purification_procedure(po, r);
C = purification_procedure2(po);
%% ����� ���������� �������
X = X_for_mixed(r, x);

%% ����� ���������
Lambda1 = Lambda_for_mixed(X);

%% �����������
lambda = zeros(m,1);
for k=1:m
    lambda(k,1) = C' * Lambda1(:,:,k) * C;
end

%% ���������� ��������� �������
% m - ����� ����� � ���������� �������, 
% N0 - ����������� ������������ ���������� �������
% s - �������� ����������� ������������
% r - ����
% X = zeros(m, r*s, r);
[m, r_s, r] = size(X);
s = fix(r_s/r);

% ��������, ����������� ��� ���������� ������� ����������
c = [real(C); imag(C)];

% %% ������������� ���������� ������� ����� ����
% phi = zeros(N,1);
% theta = zeros(N,1);
% s = zeros(N,1);
% for k=1:N
%     phi(k) = phase(X(k,2))-phase(X(k,1));
% 	theta(k) = 2*acos(abs(X(k,1)));
%     s(k) = (phase(X(k,1))+phase(X(k,2)))/2;
% end


%% ������������ ������� ����������
% ������� ����������
H = zeros(2*r_s,2*r_s);
s3=2;
% ���������� �������
X2 = zeros(m, s3, 2*r*s, r);
for j=1:m
	for l=1:r
        X2(j,:,:,l) = [real(X(j,:,l)), -imag(X(j,:,l));
                       imag(X(j,:,l)),  real(X(j,:,l))];
    end
end

% ���������
Lambda = zeros(2*r*s,2*r*s,m);

for j=1:m
	for l=1:r
        X_tmp = zeros(s,2*r*s);
        for s1=1:s3
            for s2=1:2*r*s
                X_tmp(s1,s2) = X2(j,s1,s2,l);
            end
        end
        Lambda(:,:,j) = Lambda(:,:,j) + X_tmp'*X_tmp;
    end
end

for k=1:m
   if( lambda(k) > 10^(-30))
         H = H + t(k)/lambda(k)*(Lambda(:,:,k)*c)*(Lambda(:,:,k)*c)';
    end
end

H = 2 * H;

%% ���������� �������� ����������
S = eig(H);
S = sort(S);
nu = (2*s-r)*r;
S = S(length(S)-nu+1:end-1);
d = 1./(2.*S);
d = real(d);
% d = zeros(nu, 1);
% for k=1:length(d)
%    d(k) = 1 / (2 * S(k));
% end

%% ���������
dF = sum(d);
if(isinf(dF))
    dF
end
