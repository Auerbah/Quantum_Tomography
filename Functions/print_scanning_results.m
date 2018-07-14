function print_scanning_results(filename, L_X_Cube)
maxC_Cube = max(max(L_X_Cube));
minC_Cube = min(min(L_X_Cube));
fprintf('%s\n', filename);
fprintf('Lmin = %.3f, Lmax = %.3f\n', minC_Cube, maxC_Cube);
toc
fprintf('\n');
end

