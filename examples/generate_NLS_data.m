% GENERATE_NLS_DATA  Solve the Nonlinear Schrödinger (NLS) equation
%                    using a split-step Fourier method via ode45.
%
%   Equation:  i*psi_t + (1/2)*psi_xx + |psi|^2 * psi = 0
%   Initial condition:  psi(x, 0) = 2*sech(x)   (2-soliton solution)
%
%   This script generates the spatiotemporal data (psi_sol) that is used
%   as input to the TASC-DMD benchmark in Schrodinger_equ_TASCDMD.m.
%
%   Output variable saved to workspace:
%       psi_sol  - Complex field  [nt x n_points]  |psi(x,t)|
%       x        - Spatial grid  [1 x n_points]
%       t        - Time vector   [nt x 1]
%
%   Dependencies:
%       nls_rhs.m  (must be on MATLAB path)
%
%   Reference:
%       Arani et al., "A Novel Spatiotemporal Decomposition and
%       Identification of Sparse Equations for Human Brain Deformation."

clc; clear; close all;

% -------------------------------------------------------------------------
% Spatial domain
% -------------------------------------------------------------------------
x_min    = -5;
x_max    =  5;
n_points = 500;

x_full = linspace(x_min, x_max, n_points + 1);
x      = x_full(1:n_points);          % Periodic domain (drop repeated endpoint)

% -------------------------------------------------------------------------
% Time domain
% -------------------------------------------------------------------------
t = linspace(0, pi, 21)';             % 21 snapshots over [0, pi]

% -------------------------------------------------------------------------
% Initial condition: 2-soliton  psi(x,0) = 2*sech(x)
% -------------------------------------------------------------------------
psi0 = 2 * sech(x)';                  % Column vector [n_points x 1]

% Fourier wave numbers (consistent with periodic BC on [x_min, x_max])
k = (2*pi / (x_max - x_min)) * [0 : n_points/2-1,  -n_points/2 : -1]';

% Transform initial condition to Fourier space
psit0 = fft(psi0);

% -------------------------------------------------------------------------
% Integrate in Fourier space using ode45
% -------------------------------------------------------------------------
opts = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);
[t, psit_sol] = ode45(@(tt, pp) nls_rhs(tt, pp, [], k), t, psit0, opts);

% Transform back to physical space
psi_sol = zeros(size(psit_sol));
for j = 1 : length(t)
    psi_sol(j, :) = ifft(psit_sol(j, :));
end

% -------------------------------------------------------------------------
% Visualization: |psi(x,t)|
% -------------------------------------------------------------------------
[Tgrid, Xgrid] = meshgrid(t, x);

figure;
surf(Xgrid, Tgrid, abs(psi_sol)', 'EdgeColor', 'none');
colormap('cool');
xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('$t$', 'Interpreter', 'latex', 'FontSize', 18);
zlabel('$|\psi(x,t)|$', 'Interpreter', 'latex', 'FontSize', 18);
title('Nonlinear Schrödinger Equation — Original (no noise)', 'FontSize', 14);
colorbar;
set(gca, 'FontSize', 13, 'FontWeight', 'bold');

% Save data for use in the TASC-DMD script
save('NLS_data.mat', 'psi_sol', 'x', 't', 'Xgrid', 'Tgrid');
fprintf('Data saved to NLS_data.mat\n');
