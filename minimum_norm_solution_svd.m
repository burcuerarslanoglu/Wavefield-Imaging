% Function for Part 19

function x_mn_with_svd = minimum_norm_solution_svd(A_current,chi_rho_vector)
    u_sc_current = A_current * chi_rho_vector;
    [U, S, V] = svd(A_current);
    M_current = size(A_current,1);
    U_R = U(:,1:M_current);
    S_R = S(:,1:M_current);
    V_R = V(:,1:M_current);
    Sigma_inv = diag(1 ./ diag(S_R));  %Inverse (^-1) of the diagonal elements of S
    x_mn_with_svd = V_R * Sigma_inv * U_R' * u_sc_current;
end