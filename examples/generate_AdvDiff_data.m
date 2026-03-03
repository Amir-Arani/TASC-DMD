% GENERATE_ADVDIFF_DATA  Solve the linear advection-diffusion equation
%                        using spectral discretization + ode45.
%
%   NOTE: This example is NOT included in the paper. It is provided as
%   an additional benchmark to demonstrate TASC-DMD on a PDE with
%   analytically known eigenvalues.
%
%   Equation:  u_t + U0*u_x = nu*u_xx
%
%   Analytical eigenvalues (exact, known):
%       omega_n = -nu*k_n^2  -  i*U0*k_n
%   where k_n = 2*pi*n/L
%
%   These lie on the parabola:
%       Re(omega) = -(nu/U0^2) * Im(omega)^2
%
%   Output saved to AdvDiff_data.mat:
%       u_sol        [n_points x n_snap]
%       x            [n_points x 1]
%       t            [n_snap x 1]
%       Xgrid, Tgrid
%       omega_exact  [n_points x 1]  -- analytical eigenvalues for overlay

clc; clear; close all;

% -----------------------------------------------------------------------
% Parameters
% -----------------------------------------------------------------------
U0 = 1.5;     % Advection speed
nu = 0.05;    % Viscosity
L  = 2*pi;    % Domain length

% -----------------------------------------------------------------------
% Spatial grid  (periodic BCs)
% -----------------------------------------------------------------------
n_points = 128;
x        = linspace(0, L, n_points+1)';
x        = x(1:n_points);
dx       = x(2) - x(1);

% Fourier wave numbers
n_vec = [0 : n_points/2-1,  -n_points/2 : -1]';
k     = (2*pi/L) * n_vec;

% -----------------------------------------------------------------------
% Analytical eigenvalues (ground truth for overlay)
% -----------------------------------------------------------------------
omega_exact = -nu .* k.^2  -  1i .* U0 .* k;

fprintf('=== Advection-Diffusion Benchmark (NOT in paper) ===\n');
fprintf('U0 = %.2f,  nu = %.4f\n', U0, nu);
fprintf('Analytical parabola: Re(w) = -%.4f * Im(w)^2\n\n', nu/U0^2);

% -----------------------------------------------------------------------
% Time grid  (41 snapshots)
% -----------------------------------------------------------------------
t_span = linspace(0, 4.0, 41)';

% -----------------------------------------------------------------------
% Initial condition  (multi-mode sinusoids)
% -----------------------------------------------------------------------
u0 = 0.8*sin(x) + 0.5*sin(2*x) + 0.3*sin(3*x) + 0.15*sin(4*x);

% -----------------------------------------------------------------------
% Spectral RHS  (exact linear operator in Fourier space)
% -----------------------------------------------------------------------
L_op        = -1i .* U0 .* k  -  nu .* k.^2;
advdiff_rhs = @(t, u) real(ifft(L_op .* fft(u)));

% -----------------------------------------------------------------------
% Solve with ode45
% -----------------------------------------------------------------------
fprintf('Solving with ode45...\n');
opts   = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);
[t, U] = ode45(advdiff_rhs, t_span, u0, opts);
fprintf('Done.  %d snapshots.  u range: [%.4f, %.4f]\n\n', ...
        length(t), min(U(:)), max(U(:)));

u_sol = U';     % [n_points x n_snap]

% -----------------------------------------------------------------------
% Meshgrid
% -----------------------------------------------------------------------
[Tgrid, Xgrid] = meshgrid(t, x);

% -----------------------------------------------------------------------
% Figure 1 — 3D surf  (NLS style)
% -----------------------------------------------------------------------
figure(1); clf;
surf(Xgrid, Tgrid, u_sol, 'EdgeColor', 'none');
colormap('cool');  colorbar;
caxis([min(u_sol(:)), max(u_sol(:))]);
xlabel('$x$',      'Interpreter', 'latex', 'FontSize', 18, 'FontWeight', 'bold');
ylabel('$t$',      'Interpreter', 'latex', 'FontSize', 18, 'FontWeight', 'bold');
zlabel('$u(x,t)$', 'Interpreter', 'latex', 'FontSize', 18, 'FontWeight', 'bold');
title('Advection-Diffusion — original (no noise)', 'FontSize', 14);
set(gca, 'FontSize', 13, 'FontWeight', 'bold');
view([-40 30]);
drawnow;

% -----------------------------------------------------------------------
% Figure 2 — profiles at selected times
% -----------------------------------------------------------------------
figure(2); clf;
idx  = [1, 11, 21, 31, 41];
cmap = lines(length(idx));
for k_idx = 1 : length(idx)
    plot(x, u_sol(:, idx(k_idx)), 'Color', cmap(k_idx,:), ...
         'LineWidth', 2, 'DisplayName', sprintf('t = %.2f', t(idx(k_idx))));
    hold on;
end
xlabel('$x$',      'Interpreter', 'latex', 'FontSize', 16, 'FontWeight', 'bold');
ylabel('$u(x,t)$', 'Interpreter', 'latex', 'FontSize', 16, 'FontWeight', 'bold');
title('Advection-Diffusion — traveling and decaying profiles', 'FontSize', 13);
legend('Location', 'northeast', 'FontSize', 11);
grid on;  set(gca, 'FontSize', 13, 'FontWeight', 'bold');
drawnow;

% -----------------------------------------------------------------------
% Figure 3 — analytical eigenvalue spectrum (preview)
% -----------------------------------------------------------------------
Im_range = linspace(min(imag(omega_exact)), max(imag(omega_exact)), 300);
Re_para  = -(nu/U0^2) .* Im_range.^2;

figure(3); clf;
plot(Re_para, Im_range, 'k--', 'LineWidth', 2, ...
     'DisplayName', 'Analytical parabola');
hold on;
scatter(real(omega_exact), imag(omega_exact), 40, 'k', 'filled', ...
        'DisplayName', 'Exact eigenvalues');
xlabel('Re($\omega$)', 'Interpreter', 'latex', 'FontSize', 15, 'FontWeight', 'bold');
ylabel('Im($\omega$)', 'Interpreter', 'latex', 'FontSize', 15, 'FontWeight', 'bold');
title('Analytical eigenvalue spectrum (ground truth)', 'FontSize', 13);
legend('Location', 'northeast', 'FontSize', 11);
grid on;  axis square;
set(gca, 'FontSize', 12, 'FontWeight', 'bold');
drawnow;

% -----------------------------------------------------------------------
% Save
% -----------------------------------------------------------------------
save('AdvDiff_data.mat', 'u_sol', 'x', 't', 'Xgrid', 'Tgrid', ...
     'omega_exact', 'U0', 'nu', 'L', 'k');
fprintf('Saved AdvDiff_data.mat  [%d x %d]\n', n_points, length(t));
fprintf('NOTE: This example is not included in the paper.\n');
