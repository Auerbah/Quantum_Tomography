function [filename_mat, X_matrix] = protocols_for_scanning()
filename_mat = {
%                   'X_thomson_4'
%                   'X_thomson_5'
%                   'X_thomson_6'
%                   'X_thomson_7'
%                   'X_thomson_8'
%                   'X_thomson_9'
%                   'X_thomson_10'
%                   'X_thomson_11'
%                   'X_thomson_12'
%                   'X_thomson_13'
%                   'X_thomson_14'
%                   'X_thomson_15'
%                   'X_thomson_16'
%                   'X_thomson_17'
%                   'X_thomson_18'
%                   'X_thomson_19'
%                   'X_thomson_20'
%                   'X_thomson_21'
%                   'X_thomson_22'
%                   'X_thomson_23'
%                   'X_thomson_24'
%                   'X_thomson_25'
%                   'X_thomson_26'
%                   'X_thomson_27'
%                   'X_thomson_28'
%                   'X_thomson_29'
%                   'X_thomson_30'
%                   'X_thomson_31'
%                   'X_thomson_32'
%                   'X_thomson_33'
%                   'X_thomson_34'
%                   'X_thomson_35'
%                   'X_thomson_36'
%                 'pack__1qb_4st'
%                 'pack__1qb_5st'
%                 'pack__1qb_6st'
%                 'pack__1qb_7st'
%                 'pack__1qb_8st'
%                 'pack__1qb_9st'
%                 'pack__1qb_10st'
%                 'pack__1qb_11st'
%                 'pack__1qb_12st'
%                 'pack__1qb_24st'
%                 'pack__1qb_72st'
%                 'pack__1qb_100st'
                'X_Tetrahedron'
%                 'X_Cube'
%                 'X_Octahedron'
%                 'X_Dodecahedron'
%                 'X_Icosahedron'
%                 'X_Fullerene'
%                 'X_Fullerene_dual'
                };
N_L = length(filename_mat);
X_matrix = {};
for k=1:N_L
    S = load(char(filename_mat(k)));
    X_matrix{k} = S.X;
end
end
