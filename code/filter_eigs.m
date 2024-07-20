function lambda = filter_eigs(e, omega)
    tol = 1e-2;
    condition1 = real(e) > 0 & real(e) < 1 & abs(imag(e)) < 1;
    e = e(condition1);
    
    % Sort the eigenvalues
    
    [~, condition2] = sort(imag(e), 'descend');
    e = e(condition2); 
    
    % Ensure that the last eigenvalue is the correct one
    
    condition3 = real(e) > omega + tol; 
    e = e(condition3);
    
    lambda = e(end);

end