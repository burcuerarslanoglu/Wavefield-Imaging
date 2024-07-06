close all
clear all
clc

%% Part 1

n = 0;  % zero-order 
k = 2; % second kind
z = linspace(-15,15,500);  % x-axis

H0_2 = besselh(n, k, abs(z)); % zero-order Hankel function of the second kind

% Plotting
figure
plot(z, real(H0_2), 'b-')
hold on
plot(z, imag(H0_2), 'r--')
xlabel('$z$', 'Interpreter', 'latex')
ylabel('$H_0^{(2)}(z)$', 'Interpreter', 'latex')
legend('Real Part', 'Imaginary Part')
title('Zero-Order Hankel Function of the Second Kind, $H_0^{(2)}(z)$', 'Interpreter', 'latex')
grid on
grid minor
hold off

%% Part 2 --> REPORT

%% Part 3 

k_b = 1; % for simplicity
%k_b = 2; % for a different k_b (wave number is bigger)
lambda = 2 * pi / k_b;  % wavelength of the background field

% Object domain D location (rectangle)
x_rect_obj = [0, lambda, lambda, 0, 0];
y_rect_obj = [0, 0, lambda, lambda, 0];

% Source location
rho_s = [lambda/2, 10*lambda];  
%rho_s = [lambda/2, 2*lambda]; % for a different rho_s (source is closer to the object domain)

% Plotting
figure
hold on
line(x_rect_obj, y_rect_obj, 'Color', 'b', 'LineWidth', 2) % Object domain D (rectangle) location
plot(rho_s(1), rho_s(2), 'ro', 'MarkerSize', 5, 'DisplayName', 'Source') % Source location
set(gca, 'XAxisLocation', 'top')
set(gca, 'YDir', 'reverse') % y-axis should point down
axis equal;
xlabel('x', 'Interpreter', 'latex')
xlim([0 70])
ylabel('y', 'Interpreter', 'latex')
ylim([0 70])
%title('Sketch of the Configuration')
lgd = legend('Object Domain, $\mathcal{D}$','Source, $\mathbf{\rho}_s = (\lambda/2)\mathbf{i}_x + 10\lambda\mathbf{i}_y$', 'Interpreter', 'latex');
lgd.FontSize = 14;
grid on
grid minor
%text(x_rect(2), y_rect(3), ' \leftarrow Object domain \mathcal{D}', 'Interpreter', 'latex', 'VerticalAlignment', 'bottom', 'FontSize', 13);
%text(rho_s(1), rho_s(2), ' \leftarrow Source (\lambda/2, 10\lambda)', 'Interpreter', 'latex', 'VerticalAlignment', 'bottom', 'FontSize', 13);
text(x_rect_obj(2) , y_rect_obj(3)+1, ' \leftarrow (\lambda, \lambda)', 'VerticalAlignment', 'bottom', 'FontSize', 13);
text(rho_s(1), rho_s(2)+1.5, ' \leftarrow (\lambda/2, 10\lambda)', 'VerticalAlignment', 'bottom', 'FontSize', 13);
text(10, 35, 'Wavelength of the Background Field: $\lambda = 2\pi/k_b$', 'Interpreter', 'latex', 'VerticalAlignment', 'bottom', 'FontSize', 13);
hold off

%% Part 4 --> DONE

%% Part 5

step_size = lambda/20; % Step size (uniform step size)
n_for_dim = lambda / step_size; % # of grid points along each dimension (x and y)
N = n_for_dim * n_for_dim; % Total # of grid points
disp(['Total number of grid points, N: ', num2str(N)]);

%% Part 6

x = (step_size/2):step_size:(lambda - step_size/2);
y = (step_size/2):step_size:(lambda - step_size/2);
[X, Y] = meshgrid(x, y); % the grid inside the object domain

% Computation of the incident field (eq 1 from the project manual)
difference = sqrt((X - rho_s(1)).^2 + (Y - rho_s(2)).^2); % X and Y are matrices 
u_inc = -j / 4 * besselh(n, k,  k_b * abs(difference)); % where n=0 and k=2

% Store as a 2D array
u_inc_real = real(u_inc);
u_inc_imag = imag(u_inc);
u_inc_abs = abs(u_inc);

% Plotting
figure

subplot(1, 3, 1)
imagesc([0 lambda], [0 lambda], u_inc_real)
axis equal tight
colorbar
xlabel('x', 'Interpreter', 'latex')
ylabel('y', 'Interpreter', 'latex')
title('Real part of $\hat{u}^{inc}$', 'Interpreter', 'latex')
set(gca, 'XAxisLocation', 'top')
set(gca, 'YDir', 'reverse')  % y-axis should point down

subplot(1, 3, 2)
imagesc([0 lambda], [0 lambda], u_inc_imag)
axis equal tight
colorbar
xlabel('x', 'Interpreter', 'latex')
ylabel('y', 'Interpreter', 'latex')
title('Imaginary part of $\hat{u}^{inc}$', 'Interpreter', 'latex')
set(gca, 'XAxisLocation', 'top')
set(gca, 'YDir', 'reverse')  % y-axis should point down

subplot(1, 3, 3);
imagesc([0 lambda], [0 lambda], u_inc_abs)
axis equal tight
colorbar
xlabel('x', 'Interpreter', 'latex')
ylabel('y', 'Interpreter', 'latex')
title('Absolute value of $\hat{u}^{inc}$', 'Interpreter', 'latex')
set(gca, 'XAxisLocation', 'top')
set(gca, 'YDir', 'reverse')  % y-axis should point down

%% Part 7 --> DONE ABOVE

%% Part 8 

% Arbitrary circular object definition
center_x = lambda / 2;
center_y = lambda / 2;
radius = lambda / 2;

% Initialization 
k_rho = k_b * ones(size(X));  

% Boundaries and their values
distance = sqrt((X - center_x).^2 + (Y - center_y).^2);
k_rho(distance <= radius) = 1.5 * k_b; 
k_rho(distance <= radius*0.5) = 3 * k_b;  
k_rho(distance <= radius*0.25) = 5 * k_b;  

% % Adding an asymmetric feature (dot) on the left side of the circle
% dot_center_x = center_x - radius * 0.3;  % Shift the center of the dot to the left
% dot_center_y = center_y;                % Keep the dot center aligned vertically
% dot_radius = radius * 0.1;              % Smaller radius for the dot
% 
% % Calculate distances for the dot
% dot_distance = sqrt((X - dot_center_x).^2 + (Y - dot_center_y).^2);
% k_rho(dot_distance <= dot_radius) = 10 * k_b;  % Set a high wave number for the dot

% Computingthe the contrast function, chi(rho)
chi_rho = (k_rho / k_b).^2 - 1;

%% Part 9 

% Plotting
figure
imagesc(x, y, chi_rho);
axis equal tight
colorbar
xlabel('x', 'Interpreter', 'latex');
ylabel('y', 'Interpreter', 'latex');
title('Contrast function $\chi(\rho)$', 'Interpreter', 'latex');
set(gca, 'XAxisLocation', 'top');
set(gca, 'YDir', 'reverse');  % y-axis should point down

%% Part 10

% Receiver domain
x_rect_rec_endpoints = [-lambda, 2*lambda];
%x_rect_rec_endpoints = [-2*lambda, 4*lambda];
y_rect_rec_endpoints = [1.5*lambda, 1.5*lambda];
%y_rect_rec_endpoints = [4*lambda, 4*lambda];
% Plotting
figure
hold on
line(x_rect_obj, y_rect_obj, 'Color', 'b', 'LineWidth', 2) % Object domain D (rectangle) location
plot(rho_s(1), rho_s(2), 'ro', 'MarkerSize', 5) % Source location
line(x_rect_rec_endpoints, y_rect_rec_endpoints, 'Color', 'g', 'LineWidth', 1.5) % Receiver domain D_rec location
set(gca, 'XAxisLocation', 'top')
set(gca, 'YDir', 'reverse') % y-axis should point down
axis equal
xlabel('x', 'Interpreter', 'latex')
xlim([-30 70])
ylabel('y', 'Interpreter', 'latex')
ylim([0 70])
title('Sketch of the Configuration with the Receiver Domain')
lgd = legend('Object Domain, $\mathcal{D}$','Source, $\mathbf{\rho}_s = (\lambda/2)\mathbf{i}_x + 10\lambda\mathbf{i}_y$','Receiver Domain, $\mathcal{D}^{rec}$', 'Interpreter', 'latex','Location','southeast');
lgd.FontSize = 14;
grid on
grid minor
text(x_rect_obj(2) , y_rect_obj(3)+1, ' \leftarrow (\lambda, \lambda)', 'VerticalAlignment', 'bottom', 'FontSize', 13);
text(rho_s(1), rho_s(2)+1.5, ' \leftarrow (\lambda/2, 10\lambda)', 'VerticalAlignment', 'bottom', 'FontSize', 13);
text(x_rect_rec_endpoints(2), y_rect_rec_endpoints(2)+1.5, ' \leftarrow (2\lambda, 1.5\lambda)', 'VerticalAlignment', 'bottom', 'FontSize', 13);
text(x_rect_rec_endpoints(1)-18, y_rect_rec_endpoints(1)+1.4, ' (-\lambda, 1.5\lambda) \rightarrow', 'VerticalAlignment', 'bottom', 'FontSize', 13);
text(10, 35, 'Wavelength of the Background Field: $\lambda = 2\pi/k_b$', 'Interpreter', 'latex', 'VerticalAlignment', 'bottom', 'FontSize', 13);
hold off

%% Part 11

M = 35;  % M uniformly spaced gridpoints

% Receiver domain gridpoints
x_rect_rec_gridpoints = linspace(x_rect_rec_endpoints(1), x_rect_rec_endpoints(2), M);
y_rect_rec_gridpoints = linspace(y_rect_rec_endpoints(1), y_rect_rec_endpoints(2), M);

% Plotting
figure
hold on
hold on
line(x_rect_obj, y_rect_obj, 'Color', 'b', 'LineWidth', 2) % Object domain D (rectangle) location
plot(rho_s(1), rho_s(2), 'ro', 'MarkerSize', 5) % Source location
line(x_rect_rec_endpoints, y_rect_rec_endpoints, 'Color', 'g', 'LineWidth', 1.5) % Receiver domain D_rec location
plot(x_rect_rec_gridpoints, y_rect_rec_gridpoints, 'mx', 'MarkerSize', 5) % Receiver grid points
set(gca, 'XAxisLocation', 'top')
set(gca, 'YDir', 'reverse') % y-axis should point down
axis equal;
xlabel('x', 'Interpreter', 'latex')
xlim([-30 70])
ylabel('y', 'Interpreter', 'latex')
ylim([0 70])
%title('Sketch of the Configuration with the Grid Points of the Receiver Domain')
lgd = legend('Object Domain, $\mathcal{D}$','Source, $\mathbf{\rho}_s = (\lambda/2)\mathbf{i}_x + 10\lambda\mathbf{i}_y$','Receiver Domain, $\mathcal{D}^{rec}$','Receiver Grid Points, $M=35$', 'Interpreter', 'latex','Location','southeast');
lgd.FontSize = 14;
grid on
grid minor
text(x_rect_obj(2) , y_rect_obj(3)+1, ' \leftarrow (\lambda, \lambda)', 'VerticalAlignment', 'bottom', 'FontSize', 13);
text(rho_s(1), rho_s(2)+1.5, ' \leftarrow (\lambda/2, 10\lambda)', 'VerticalAlignment', 'bottom', 'FontSize', 13);
text(x_rect_rec_endpoints(2), y_rect_rec_endpoints(2)+1.5, ' \leftarrow (2\lambda, 1.5\lambda)', 'VerticalAlignment', 'bottom', 'FontSize', 13);
text(x_rect_rec_endpoints(1)-18, y_rect_rec_endpoints(1)+1.4, ' (-\lambda, 1.5\lambda) \rightarrow', 'VerticalAlignment', 'bottom', 'FontSize', 13);
text(10, 35, 'Wavelength of the Background Field: $\lambda = 2\pi/k_b$', 'Interpreter', 'latex', 'VerticalAlignment', 'bottom', 'FontSize', 13);
hold off

% Plotting (for 20)
% figure
% hold on
% hold on
% line(x_rect_obj, y_rect_obj, 'Color', 'b', 'LineWidth', 2) % Object domain D (rectangle) location
% plot(rho_s(1), rho_s(2), 'ro', 'MarkerSize', 5) % Source location
% line(x_rect_rec_endpoints, y_rect_rec_endpoints, 'Color', 'g', 'LineWidth', 1.5) % Receiver domain D_rec location
% plot(x_rect_rec_gridpoints, y_rect_rec_gridpoints, 'mx', 'MarkerSize', 5) % Receiver grid points
% set(gca, 'XAxisLocation', 'top')
% set(gca, 'YDir', 'reverse') % y-axis should point down
% axis equal;
% xlabel('x', 'Interpreter', 'latex')
% xlim([-200 300])
% ylabel('y', 'Interpreter', 'latex')
% ylim([0 300])
% %title('Sketch of the Configuration with the Grid Points of the Receiver Domain')
% lgd = legend('Object Domain, $\mathcal{D}$','Source, $\mathbf{\rho}_s = (\lambda/2)\mathbf{i}_x + 10\lambda\mathbf{i}_y$','Receiver Domain, $\mathcal{D}^{rec}$','Receiver Grid Points, $M=10$', 'Interpreter', 'latex','Location','southeast');
% lgd.FontSize = 14;
% grid on
% grid minor
% text(x_rect_obj(2) , y_rect_obj(3)+1.2, ' \leftarrow (\lambda, \lambda)', 'VerticalAlignment', 'bottom', 'FontSize', 13);
% text(rho_s(1), rho_s(2)+1.5, ' \leftarrow (\lambda/2, 10\lambda)', 'VerticalAlignment', 'bottom', 'FontSize', 13);
% text(x_rect_rec_endpoints(2), y_rect_rec_endpoints(2)+1.3, ' \leftarrow (4\lambda, 4\lambda)', 'VerticalAlignment', 'bottom', 'FontSize', 13);
% text(x_rect_rec_endpoints(1)-65, y_rect_rec_endpoints(1)+1.3, ' (-2\lambda, 4\lambda) \rightarrow', 'VerticalAlignment', 'bottom', 'FontSize', 13);
% text(10, 55, 'Wavelength of the Background Field: $\lambda = 2\pi/k_b$', 'Interpreter', 'latex', 'VerticalAlignment', 'bottom', 'FontSize', 13);
% hold off

%% Part 12 --> REPORT

%% Part 13

number_elements = numel(chi_rho);  % Total number of elements in the 2D array of the contrast function
chi_rho_vector = reshape(chi_rho, [number_elements, 1]); % Reshaping of chi_rho in one big column vector

%% Part 14 

A_initial_parameters = build_system_matrix(X, Y, chi_rho_vector, u_inc, x_rect_rec_endpoints, y_rect_rec_endpoints, n, k, k_b, M);

%% Part 15 --> Emin değilim.
 
% Number of receiver locations (vary M to investigate its effect)
M_values = [10, 20, 30, 40, 50];  % Different values of M to investigate
singular_values = cell(length(M_values), 1);

% Preallocate A for different values of M
A_all = cell(length(M_values), 1);

for i = 1:length(M_values)
    M_current = M_values(i);
    % Build the system matrix A for current M
    A = build_system_matrix(X, Y, chi_rho_vector, u_inc, x_rect_rec_endpoints, y_rect_rec_endpoints, n, k, k_b, M_current);
    
    % Store the system matrix
    A_all{i} = A;
    
    % Compute the singular value decomposition
    [U, S, V] = svd(A);
    
    % Store the singular values
    singular_values{i} = diag(S);
end

% Plot the singular values for different M
figure
hold on;
num_rows = 2;
num_cols = length(M_values);
for i = 1:length(M_values)
    % Regular
    subplot(num_rows, num_cols, i);
    plot(singular_values{i}, '-', 'Linewidth', 1.5, 'DisplayName', sprintf('M = %d', M_values(i)));
    title(sprintf('M = %d', M_values(i)));
    xlabel('Number of Singular Values');
    ylabel('Singular Values');
    xlim([1  M_values(i)])
    grid on;
  
    % 10log
    subplot(num_rows, num_cols, i + num_cols);
    plot(log10(singular_values{i}), '-', 'Linewidth', 1.5, 'DisplayName', sprintf('M = %d', M_values(i)));
    title(sprintf('M = %d', M_values(i)));
    xlabel('Number of Singular Values');
    ylabel('10log of Singular values');
    xlim([1  M_values(i)])
    grid on;
end

%% Part 16

% Scattered field with the initial parameters (M=30 in Part 11)
u_sc_initial_parameters = A_initial_parameters * chi_rho_vector;

%% Part 17

% Minimum norm solution using SVD (eq. 4.29)
[U, S, V] = svd(A_initial_parameters);
U_R = U(:,1:M);
S_R = S(:,1:M);
V_R = V(:,1:M);
Sigma_inv = diag(1 ./ diag(S_R));  %Inverse (^-1) of the diagonal elements of S
x_mn_svd = V_R * Sigma_inv * U_R' * u_sc_initial_parameters;

% Minimum norm solution using pinv (eq. 4.30)
x_mn_pinv = pinv(A_initial_parameters) * u_sc_initial_parameters;

%% Part 18 
x_mn_2D_with_svd = reshape(x_mn_svd, [size(chi_rho,1), size(chi_rho,2)]);
x_mn_2D_with_pinv = reshape(x_mn_pinv, [size(chi_rho,1), size(chi_rho,2)]);

figure
sgtitle('Comparison Plots')

subplot(1, 3, 1);
imagesc(x,y,chi_rho);
colorbar; 
title('Original $\chi(\rho)$', 'Interpreter', 'latex');
axis equal tight
xlabel('x', 'Interpreter', 'latex');
ylabel('y', 'Interpreter', 'latex');
set(gca, 'XAxisLocation', 'top');
set(gca, 'YDir', 'reverse');  % y-axis should point down

subplot(1, 3, 2);
imagesc(x,y,real(x_mn_2D_with_svd));
colorbar; 
title('Reconstructed $\hat{\chi}(\rho)$ with SVD', 'Interpreter', 'latex');
axis equal tight
xlabel('x', 'Interpreter', 'latex');
ylabel('y', 'Interpreter', 'latex');
set(gca, 'XAxisLocation', 'top');
set(gca, 'YDir', 'reverse');  % y-axis should point down

subplot(1, 3, 3);  
imagesc(x,y,real(x_mn_2D_with_pinv));
colorbar;  
title('Reconstructed $\hat{\chi}(\rho)$ with pinv(.)', 'Interpreter', 'latex');
axis equal tight 
xlabel('x', 'Interpreter', 'latex');
ylabel('y', 'Interpreter', 'latex');
set(gca, 'XAxisLocation', 'top');
set(gca, 'YDir', 'reverse');  % y-axis should point down

figure
sgtitle('Reconstructed $\hat{\chi}(\rho)$ with SVD imaginary part and absolute value.', 'Interpreter', 'latex');

subplot(1, 2, 1);
imagesc(x,y,imag(x_mn_2D_with_svd));
colorbar; 
title('Imaginary Part of the Reconstructed $\hat{\chi}(\rho)$ with SVD', 'Interpreter', 'latex');
axis equal tight
xlabel('x', 'Interpreter', 'latex');
ylabel('y', 'Interpreter', 'latex');
set(gca, 'XAxisLocation', 'top');
set(gca, 'YDir', 'reverse');  % y-axis should point down

subplot(1, 2, 2);
imagesc(x,y,abs(x_mn_2D_with_svd));
colorbar; 
title('Absolute Value of the Reconstructed $\hat{\chi}(\rho)$ with SVD', 'Interpreter', 'latex');
axis equal tight
xlabel('x', 'Interpreter', 'latex');
ylabel('y', 'Interpreter', 'latex');
set(gca, 'XAxisLocation', 'top');
set(gca, 'YDir', 'reverse');  % y-axis should point down

error_svd_mse = immse(real(x_mn_2D_with_svd), chi_rho);
error_svd_mae = mae(real(x_mn_2D_with_svd) - chi_rho);

error_pinv_mse = immse(real(x_mn_2D_with_pinv), chi_rho);
error_pinv_mae = mae(real(x_mn_2D_with_pinv) - chi_rho);

error_table = table({'MSE'; 'MAE'}, [error_svd_mse; error_svd_mae], [error_pinv_mse; error_pinv_mae], ...
    'VariableNames', {'ErrorType', 'SVD', 'PINV'});

disp(error_table);

%% Part 19 (for different numbers of M, produce images, in subplots)

M_values = [10, 20, 30, 40, 50];

% Initialization
x_mn_2D_with_svd = cell(length(M_values), 1);

for i = 1:length(M_values)
    M_curent = M_values(i);
    % Build the system matrix A for current M
    A_current = build_system_matrix(X, Y, chi_rho_vector, u_inc, x_rect_rec_endpoints, y_rect_rec_endpoints, n, k, k_b, M_curent);

    % Minimum norm solution construction for current A
    x_mn_with_svd = minimum_norm_solution_svd(A_current,chi_rho_vector);
    
    x_mn_2D_with_svd{i} = reshape(x_mn_with_svd, [size(chi_rho,1), size(chi_rho,2)]);
    
end

% Plotting reconstructions for different M
figure
sgtitle('Reconstructed $\hat{\chi}(\rho)$s with SVD for Different Number of Receivers $M$ ', 'Interpreter', 'latex');
hold on;
for i = 1:length(M_values)
    subplot(1,length(M_values), i);
    imagesc(x,y,real(x_mn_2D_with_svd{i}));
    colorbar; 
    title(['$M = ' num2str(M_values(i)) '$'], 'Interpreter', 'latex');  % Dynamic title showing M value
    axis equal tight 
    xlabel('x', 'Interpreter', 'latex');
    ylabel('y', 'Interpreter', 'latex');
    set(gca, 'XAxisLocation', 'top');
    set(gca, 'YDir', 'reverse');  % y-axis should point down
end
hold off;

error_svd_mse = zeros(1, length(M_values));
error_svd_mae = zeros(1, length(M_values));

% Calculate errors for each M
for i = 1:length(M_values)
    error_svd_mse(i) = immse(real(x_mn_2D_with_svd{i}), chi_rho);
    error_svd_mae(i) = mae(real(x_mn_2D_with_svd{i}) - chi_rho);
end

error_table = table(M_values', error_svd_mse', error_svd_mae', ...
    'VariableNames', {'M', 'MSE_SVD', 'MAE_SVD'});

disp(error_table);

%% Part 20

% Add complex Gaussian noise
%avg_magnitude_real_part = mean(real(u_sc_initial_parameters));
%avg_magnitude_imag_part = mean(imag(u_sc_initial_parameters));
%noise_level = 0.000000000005;  % This should be adjusted based on the magnitude of u_sc
noise_level = 0.001;
noise_level = 1;
noise_real = noise_level * randn(size(u_sc_initial_parameters));
noise_imag = noise_level * randn(size(u_sc_initial_parameters));

complex_noise = noise_real + 1j * noise_imag;
u_sc_noisy = u_sc_initial_parameters + complex_noise;

% Use the noisy u_sc for image reconstruction
x_mn_pinv_noisy = pinv(A_initial_parameters) * u_sc_noisy; 

% Minimum norm solution using SVD (eq. 4.29)
[U, S, V] = svd(A_initial_parameters);
U_R = U(:,1:M);
S_R = S(:,1:M);
V_R = V(:,1:M);
Sigma_inv = diag(1 ./ diag(S_R));  %Inverse (^-1) of the diagonal elements of S
x_mn_svd_noisy = V_R * Sigma_inv * U_R' * u_sc_noisy;

% Reshaping
x_mn_2D_pinv_noisy = reshape(x_mn_pinv_noisy, [size(chi_rho,1), size(chi_rho,2)]);
x_mn_2D_svd_noisy = reshape(x_mn_svd_noisy, [size(chi_rho,1), size(chi_rho,2)]);

figure
imagesc(x,y,abs(x_mn_2D_svd_noisy));
colorbar;  
%title('Reconstructed Contrast')
axis equal tight 
xlabel('x', 'Interpreter', 'latex');
ylabel('y', 'Interpreter', 'latex');
set(gca, 'XAxisLocation', 'top');
set(gca, 'YDir', 'reverse');  % y-axis should point down

% Check condition number
condA = cond(A_initial_parameters);
fprintf('Condition Number of A: %e\n', condA);

%Tikhonov regularization
lambda = 0.01;  % Regularization parameter
A_reg = A_initial_parameters' * A_initial_parameters + lambda * eye(size(A_initial_parameters, 2));
b_reg = A_initial_parameters' * u_sc_noisy;
x_reg = A_reg \ b_reg;

x_mn_2D_reg = reshape(x_reg, size(chi_rho));

figure
imagesc(x,y,real(x_mn_2D_reg));
colorbar;
axis equal tight
%title('Regularized Reconstruction');
xlabel('x', 'Interpreter', 'latex');
ylabel('y', 'Interpreter', 'latex');
set(gca, 'XAxisLocation', 'top');
set(gca, 'YDir', 'reverse');  % y-axis should point down

