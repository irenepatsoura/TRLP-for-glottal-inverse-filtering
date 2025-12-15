function a = mytrlp(frame, p, a_prev) 
    lambda1 = 1; % Initial Regularization constant
    lambda2 = 0.9; % Initial Regularization constant
    % Autocorrelation method
    R = xcorr(frame, p);
    R = R(p+1:end);
    R_matrix = toeplitz(R(1:p));
    r_vector = R(2:p+1);
    
    % Regularization term
    a_reg = lambda2 * a_prev;
    R_reg = R_matrix + lambda1 * eye(p);
    r_reg = r_vector + lambda1 * a_reg;
    
    % Solve for new LPC coefficients
    a = R_reg \ r_reg;
end