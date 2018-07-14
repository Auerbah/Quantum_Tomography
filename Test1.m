clc
clear
close
po = [[0 0 1 0]
      [1 0 0 0]
      [0 0 0 0]
      [0 0 0 0]];
%     trace(C1*C1)
    %trace_A 
    s = 2;
    r = 2;
    po_a = zeros(s,s);
%     po_pur = C1*C1;
    for i=1:s
        po_a = po_a + kron(basis(r,i)',eye(s))*po*kron(basis(r,i),eye(s));
   end
    po_a
    
   po_b = zeros(s,s);
%  po_pur = C1*C1;
    for i=1:s
        po_b = po_b + kron(eye(s),basis(r,i)')*po*kron(eye(s),basis(r,i));
    end
   po_b