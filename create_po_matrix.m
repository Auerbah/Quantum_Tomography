function po = create_po_matrix(r,n) 

%% ќртогональные состо€ни€, на основе которых строитс€ матрица плотности
psi0 = [1; 0];
psi1 = [0; 1];

phi0 = 1 / sqrt(2) * [1; 1];
phi1 = 1 / sqrt(2) * [1; -1];

chi0 = 1 / sqrt(2) * [1; 1i];
chi1 = 1 / sqrt(2) * [1; -1i];

%% ћатрица плотности
% if r == 1
%     po = kron_S(chi0,n)*kron_S(chi0,n)';
% elseif r == n^2
%     po = kron_S(psi0,n)*kron_S(psi0,n)'*p + ...
%          kron_S(chi0,n)*kron_S(chi0,n)'*(1-p);
% else
%     po = kron_S(psi0,n)*kron_S(psi0,n)'*p + ...
%          kron_S(chi0,n)*kron_S(chi0,n)'*(1-p);
% end

po = [eye(r), zeros(r,2^n-r); zeros(2^n-r,2^n)];
po = po/trace(po);

if ~check_po_matrix(po)
	error('')
end

end

