% SCHRODINGER_EQU_TASCDMD  Benchmark TASC-DMD, Exact DMD, and fbDMD on
%                           the Nonlinear Schrödinger (NLS) equation.
%
%   This script reproduces Figures 2a–2i of:
%       [TASC DMD Paper — citation to be added after publication]

%
%   Workflow:
%       1. Load NLS data (run generate_NLS_data.m first, or load NLS_data.mat)
%       2. Optionally add noise to the clean field
%       3. Run TASC-DMD (time-augmented, space-contracted DMD)
%       4. Run Exact DMD and fbDMD on the same data
%       5. Compare eigenvalue spectra and field reconstructions
%
%   Dependencies (all in /src):
%       DMD_exact.m, fbDMD.m, nls_rhs.m
%
%   Usage:
%       1. Run generate_NLS_data.m once to produce NLS_data.mat
%       2. Run this script
%
%   Parameters you may want to tune:
%       noise_level  - Fractional noise amplitude (0 = no noise, 0.2 = 20%)
%       d            - Number of time-delay embeddings
%       r0           - Rank of the spatial basis per delay level
%       r1           - Rank used for the final DMD step

clc; clear; close all;

% -------------------------------------------------------------------------
% 0. Setup
% -------------------------------------------------------------------------
%addpath('../src');          % Add DMD function files to path

% -------------------------------------------------------------------------
% 1. Load NLS data
% -------------------------------------------------------------------------
% Run generate_NLS_data.m first if NLS_data.mat does not exist.
if ~isfile('NLS_data.mat')
    error('NLS_data.mat not found. Please run generate_NLS_data.m first.');
end
load('NLS_data.mat', 'psi_sol', 'x', 't', 'Xgrid', 'Tgrid');

dt = t(2) - t(1);

% -------------------------------------------------------------------------
% 2. Add noise (set noise_level = 0 for clean data)
% -------------------------------------------------------------------------
noise_level = 0.20;                       % 20% fractional Gaussian noise
u_clean     = abs(psi_sol)';              % [n_points x nt], real amplitude
noise_mask  = noise_level .* u_clean;
u_noisy     = u_clean + noise_mask .* randn(size(u_clean));

% Choose which dataset to analyse (swap comments to toggle noise)
X_input = u_noisy;    % noisy
%X_input = u_clean;  % clean

% -------------------------------------------------------------------------
% 3. TASC-DMD
% -------------------------------------------------------------------------
% --- Hyperparameters ---
tp1 = 1;    % First snapshot index
tp2 = 21;   % Last  snapshot index
d   = 11;   % Number of time-delay levels
r0  = 11;   % Rank of spatial basis per delay (space-contraction rank)
r1  = 9;    % Rank for the final DMD step (mode truncation)

X_win  = X_input(:, tp1:tp2);            % Windowed data [n x m]
n      = size(X_win, 1);
m      = size(X_win, 2);
dt1    = dt;

% --- Build 3D time-delay Hankel tensor and contracted Hankel matrix ---
H_tde = zeros(n, m-d+1, d);   % [n x (m-d+1) x d]
Ur0   = zeros(n, r0,    d);   % Spatial basis for each delay level
h     = [];                    % Reduced (contracted) Hankel matrix

for j = 1 : d
    H_tde(:, :, j) = X_win(:, j : m-d+j);        % Delay-j snapshot window
    [U_tmp, ~, ~]  = svd(H_tde(:, :, j), 'econ');
    Ur0(:, :, j)   = U_tmp(:, 1:r0);              % Keep leading r0 modes
    h_tde = Ur0(:, :, j)' * H_tde(:, :, j);       % Project to low-dim space
    h = [h; h_tde]; %#ok<AGROW>                   % Stack: [(r0*d) x (m-d+1)]
end

nt2 = size(h, 2);          % = m - d + 1
X1h = h(:, 1:end-1);
X2h = h(:, 2:end);

% --- DMD on contracted Hankel matrix ---
[U_dmd, S_dmd, V_dmd] = svd(X1h, 'econ');
Ur = U_dmd(:, 1:r1);
Sr = S_dmd(1:r1, 1:r1);
Vr = V_dmd(:, 1:r1);

Atilde = Ur' * X2h * Vr / Sr;
[W, D]  = eig(Atilde);

Phi    = X2h * (Vr / Sr) * W;   % TASC-DMD modes
Lambda = diag(D);
omega  = log(Lambda) / dt1;

b = Phi \ X1h(:, 1);

time_dynamics = zeros(r1, nt2);
for i = 1 : nt2
    time_dynamics(:, i) = b .* exp(omega * (i-1) * dt1);
end

% --- Map modes back to physical space ---
H_dmd = zeros(n, nt2, d);
for j = 1 : d
    Phi_rom_j      = Phi(1+(j-1)*r0 : j*r0, :);
    psi_phys_j     = Ur0(:, :, j) * Phi_rom_j;
    H_dmd(:, :, j) = psi_phys_j * time_dynamics;
end

% Concatenate delay levels into a single reconstruction matrix
X_TASC = H_dmd(:, :, 1);
for j = 2 : d
    X_TASC = [X_TASC, H_dmd(:, end, j)]; %#ok<AGROW>
end

% --- Correlation coefficient Cv (per snapshot, per delay level) ---
Cv       = zeros(nt2, d);
Cv_cmplx = zeros(nt2, d);
for j = 1 : d
    for i = 1 : nt2
        Xn_org         = normalize(H_tde(:, i, j), 'norm');
        Xn_dmd         = normalize(H_dmd(:, i, j), 'norm');
        Cv_cmplx(i, j) = dot(Xn_org, Xn_dmd);
        Cv(i, j)       = abs(Cv_cmplx(i, j));
    end
end

% -------------------------------------------------------------------------
% 4. Exact DMD (baseline)
% -------------------------------------------------------------------------
[~, ~, ~, omega_edmd, ~, X_edmd, Cv_edmd] = ...
    DMD_exact(X_input(:, tp1:tp2), dt1, 1, 1, r1, 0);

% -------------------------------------------------------------------------
% 5. Forward-Backward DMD (baseline)
% -------------------------------------------------------------------------
[~, ~, ~, omega_fbdmd, ~, X_fbdmd, Cv_fbdmd] = ...
    fbDMD(X_input(:, tp1:tp2), dt1, 1, 1, r1);

% -------------------------------------------------------------------------
% 6. Eigenvalue spectrum comparison
% -------------------------------------------------------------------------
figure;
hold on; grid on; axis square;
scatter(real(omega_edmd),  imag(omega_edmd),  50, [0.2 0.4 0.8], 'o', 'LineWidth', 2, 'DisplayName', 'Exact DMD');
scatter(real(omega_fbdmd), imag(omega_fbdmd), 50, [1.0 0.5 0.0], 's', 'LineWidth', 2, 'DisplayName', 'fbDMD');
scatter(real(omega),       imag(omega),       50, [0.8 0.1 0.1], 'd', 'LineWidth', 2, 'DisplayName', 'TASC-DMD');
xlabel('Re($\omega_k$)', 'Interpreter', 'latex', 'FontSize', 16, 'FontWeight', 'bold');
ylabel('Im($\omega_k$)', 'Interpreter', 'latex', 'FontSize', 16, 'FontWeight', 'bold');
legend('Location', 'northeastoutside');
xlim([-14  1]);  ylim([-25 25]);
xticks(-14:2:1); yticks(-25:5:25);
set(gca, 'FontSize', 13, 'FontWeight', 'bold');
title(sprintf('Eigenvalue spectrum  (r = %d, noise = %.0f%%)', r1, noise_level*100), 'FontSize', 13);

% -------------------------------------------------------------------------
% 7. Field reconstruction comparison
% -------------------------------------------------------------------------
plot_params = {'EdgeColor', 'none'};
cax_lim     = [0, 3.8];

datasets = {u_clean(:, tp1:tp2), 'Original (clean)';
            u_noisy,             'Original (noisy)';
            real(X_edmd),        'Exact DMD';
            real(X_fbdmd),       'fbDMD';
            real(X_TASC),        'TASC-DMD'};

for k = 1 : size(datasets, 1)
    Zdata = datasets{k, 1};
    ttl   = datasets{k, 2};
    if size(Zdata, 2) ~= size(Tgrid, 2)
        continue;   % Skip if time dimension does not match the grid
    end
    figure;
    surf(Xgrid, Tgrid, Zdata, plot_params{:});
    colormap('cool');  caxis(cax_lim);  colorbar;  zlim([0 4]);
    xlabel('$x$',           'Interpreter', 'latex', 'FontSize', 18, 'FontWeight', 'bold');
    ylabel('$t$',           'Interpreter', 'latex', 'FontSize', 18, 'FontWeight', 'bold');
    zlabel('$|\psi(x,t)|$', 'Interpreter', 'latex', 'FontSize', 18, 'FontWeight', 'bold');
    title(['NLS: ', ttl], 'FontSize', 14);
    set(gca, 'FontSize', 13, 'FontWeight', 'bold');
end
