% Function for Parts 14, 15, and 19 (we can change M!)

function A = build_system_matrix(X, Y, chi_rho_vector, u_inc, x_rect_rec_endpoints, y_rect_rec_endpoints, n, k, k_b, M)

    N = numel(X); % Number of gridpoints 
    A = zeros(M, N); % Initialize the system matrix A

    % Receiver domain gridpoints determined according to M
    x_rect_rec_gridpoints = linspace(x_rect_rec_endpoints(1), x_rect_rec_endpoints(2), M);
    y_rect_rec_gridpoints = linspace(y_rect_rec_endpoints(1), y_rect_rec_endpoints(2), M);

    %Reshaping everything to column vector form for compatibility (same as Part 13)
    number_elements_X = numel(X);  
    X_vector = reshape(X, [number_elements_X, 1]); 
    number_elements_Y = numel(Y);  
    Y_vector = reshape(Y, [number_elements_Y, 1]); 
    number_elements_u_inc = numel(u_inc);  
    u_inc_vector = reshape(u_inc, [number_elements_u_inc, 1]); 

    for m = 1:M % For each receiver point 
        % Coordinates of the m-th receiver
        x_current = x_rect_rec_gridpoints(m);
        y_current = y_rect_rec_gridpoints(m);

        % Distance from each point (vector) in the domain to the current receiver
        distance = sqrt((X_vector - x_current).^2 + (Y_vector - y_current).^2);

        % Green's function 
        G_hat = -j / 4 * besselh(n, k, k_b * abs(distance));

        % Scattered field computation (eq. 3 from project manual)
        A(m, :) = -(k_b^2) * (G_hat .* chi_rho_vector .* u_inc_vector);
    end
end