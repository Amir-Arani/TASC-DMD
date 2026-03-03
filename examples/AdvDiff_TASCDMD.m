% ADVDIFF_TASCDMD  Benchmark TASC-DMD, Exact DMD, and fbDMD on the
%                  linear advection-diffusion equation with KNOWN eigenvalues.
%
%   NOTE: This example is NOT included in the paper. It is provided as
%   an additional benchmark with analytically known eigenvalues, allowing
%   direct quantitative validation of each DMD method against ground truth.
%
%   Equation:  u_t + U0*u_x = nu*u_xx
%
%   Analytical eigenvalues:  omega_n = -nu*k_n^2 - i*U0*k_n
%   In the complex plane these lie on the parabola:
%       Re(omega) = -(nu/U0^2) * Im(omega)^2
%
%   The analytical parabola is overlaid on the eigenvalue plots so
%   accuracy of each method can be assessed directly against truth.
%   Tested at 50% noise. Cv reported for all three methods.
%
%   Run generate_AdvDiff_data.m first to produce AdvDiff_data.mat.
%   Dependencies:  DMD_exact.m,  fbDMD.m

clc; clear; close all;

% -----------------------------------------------------------------------
% 0. Load data
% -----------------------------------------------------------------------
if ~isfile('AdvDiff_data.mat')
    error('AdvDiff_data.mat not found. Run generate_AdvDiff_data.m first.');
end
load('AdvDiff_data.mat', 'u_sol', 'x', 't', 'Xgrid', 'Tgrid', ...
     'omega_exact', 'U0', 'nu');

dt      = t(2) - t(1);
u_clean = u_sol;

fprintf('=== Advection-Diffusion Benchmark (NOT in paper) ===\n');
fprintf('U0 = %.2f,  nu = %.4f\n', U0, nu);
fprintf('Analytical parabola: Re(w) = -%.4f * Im(w)^2\n\n', nu/U0^2);

% -----------------------------------------------------------------------
% 1. Add noise  (50%)
% -----------------------------------------------------------------------
noise_level = 0.50;
rng(42);
u_noisy = u_clean + noise_level .* abs(u_clean) .* randn(size(u_clean));

% -----------------------------------------------------------------------
% 2. TASC-DMD hyperparameters
% -----------------------------------------------------------------------
tp1 = 1;
tp2 = size(u_clean, 2);
d   = 11;
r0  = 11;
r1  = 9;

% -----------------------------------------------------------------------
% Helper: run TASC-DMD, return omega, reconstruction, and mean Cv
% -----------------------------------------------------------------------
function [omega, X_TASC, Cv_mean] = run_TASC(X, dt1, d, r0, r1)
    [n, m] = size(X);
    H_tde  = zeros(n, m-d+1, d);
    Ur0    = zeros(n, r0, d);
    h      = [];

    for j = 1:d
        H_tde(:,:,j) = X(:, j:m-d+j);
        [Utmp,~,~]   = svd(H_tde(:,:,j), 'econ');
        Ur0(:,:,j)   = Utmp(:, 1:r0);
        h = [h; Ur0(:,:,j)' * H_tde(:,:,j)]; %#ok<AGROW>
    end

    nt2 = size(h, 2);
    X1h = h(:, 1:end-1);
    X2h = h(:, 2:end);

    [U,S,V] = svd(X1h, 'econ');
    Ur = U(:,1:r1);  Sr = S(1:r1,1:r1);  Vr = V(:,1:r1);
    [W,D]   = eig(Ur'*X2h*Vr/Sr);

    Phi   = X2h*(Vr/Sr)*W;
    omega = log(diag(D)) / dt1;
    b     = Phi \ X1h(:,1);

    T_dyn = zeros(r1, nt2);
    for i = 1:nt2
        T_dyn(:,i) = b .* exp(omega*(i-1)*dt1);
    end

    H_dmd = zeros(n, nt2, d);
    for j = 1:d
        Prj          = Phi(1+(j-1)*r0 : j*r0, :);
        H_dmd(:,:,j) = (Ur0(:,:,j) * Prj) * T_dyn;
    end

    X_TASC = H_dmd(:,:,1);
    for j = 2:d
        X_TASC = [X_TASC, H_dmd(:,end,j)]; %#ok<AGROW>
    end

    Cv = zeros(nt2, d);
    for j = 1:d
        for i = 1:nt2
            a       = normalize(H_tde(:,i,j), 'norm');
            b2      = normalize(H_dmd(:,i,j), 'norm');
            Cv(i,j) = abs(dot(a, b2));
        end
    end
    Cv_mean = mean(Cv(:));
end

% -----------------------------------------------------------------------
% 3. Run all three methods on clean and noisy data
% -----------------------------------------------------------------------
X_in   = {u_clean, u_noisy};
labels = {'Clean', '50% noise'};

omega_tasc  = cell(2,1);   Xrec_tasc  = cell(2,1);
omega_edmd  = cell(2,1);   Xrec_edmd  = cell(2,1);
omega_fbdmd = cell(2,1);   Xrec_fbdmd = cell(2,1);
Cv_table    = zeros(2,3);

fprintf('%-14s  %-12s  %-12s  %-12s\n', 'Dataset', 'TASC-DMD', 'Exact DMD', 'fbDMD');
fprintf('%s\n', repmat('-',1,54));

for di = 1:2
    X = X_in{di}(:, tp1:tp2);

    [omega_tasc{di},  Xrec_tasc{di},  Cv_table(di,1)] = run_TASC(X, dt, d, r0, r1);

    [~,~,~, omega_edmd{di},  ~, Xrec_edmd{di},  Cv_vec] = DMD_exact(X, dt, 1, 1, r1, 0);
    Cv_table(di,2) = mean(Cv_vec);

    [~,~,~, omega_fbdmd{di}, ~, Xrec_fbdmd{di}, Cv_vec] = fbDMD(X, dt, 1, 1, r1);
    Cv_table(di,3) = mean(Cv_vec);

    fprintf('%-14s  %-12.4f  %-12.4f  %-12.4f\n', ...
        labels{di}, Cv_table(di,1), Cv_table(di,2), Cv_table(di,3));
end
fprintf('\nCv = mean global correlation coefficient (1 = perfect)\n\n');

% -----------------------------------------------------------------------
% Analytical parabola for overlay
% -----------------------------------------------------------------------
Im_range = linspace(-20, 20, 500);
Re_para  = -(nu/U0^2) .* Im_range.^2;

% -----------------------------------------------------------------------
% Shared surf helper
% -----------------------------------------------------------------------
function do_surf(Xg, Tg, Z, ttl, cax, cmap)
    [~, nt] = size(Tg);
    Z = real(Z);
    if size(Z,2) < nt
        Zpad = NaN(size(Xg));  Zpad(:,1:size(Z,2)) = Z;  Z = Zpad;
    elseif size(Z,2) > nt
        Z = Z(:,1:nt);
    end
    surf(Xg, Tg, Z, 'EdgeColor', 'none');
    colormap(cmap);  clim(cax);  colorbar;
    xlabel('$x$',      'Interpreter','latex','FontSize',18,'FontWeight','bold');
    ylabel('$t$',      'Interpreter','latex','FontSize',18,'FontWeight','bold');
    zlabel('$u(x,t)$', 'Interpreter','latex','FontSize',18,'FontWeight','bold');
    title(ttl, 'FontSize', 13);
    set(gca, 'FontSize', 13, 'FontWeight', 'bold');
    view([-40 30]);
end

cax  = [min(u_clean(:))*1.1,  max(u_clean(:))*1.1];
cmap = 'cool';

% -----------------------------------------------------------------------
% Figures 1-2 — original fields
% -----------------------------------------------------------------------
figure(1); clf;
do_surf(Xgrid, Tgrid, u_clean, ...
    'Advection-Diffusion — Original (clean)', cax, cmap);
drawnow;

figure(2); clf;
do_surf(Xgrid, Tgrid, u_noisy, ...
    'Advection-Diffusion — Original (50% noise)', cax, cmap);
drawnow;

% -----------------------------------------------------------------------
% Figures 3-5 — reconstructions (clean)
% -----------------------------------------------------------------------
titles_clean = {'Exact DMD — clean', 'fbDMD — clean', 'TASC-DMD — clean'};
data_clean   = {Xrec_edmd{1}, Xrec_fbdmd{1}, Xrec_tasc{1}};
for k = 1:3
    figure(2+k); clf;
    do_surf(Xgrid, Tgrid, data_clean{k}, ...
        sprintf('%s  (Cv = %.4f)', titles_clean{k}, Cv_table(1,k)), cax, cmap);
    drawnow;
end

% -----------------------------------------------------------------------
% Figures 6-8 — reconstructions (noisy)
% -----------------------------------------------------------------------
titles_noisy = {'Exact DMD — 50% noise', 'fbDMD — 50% noise', 'TASC-DMD — 50% noise'};
data_noisy   = {Xrec_edmd{2}, Xrec_fbdmd{2}, Xrec_tasc{2}};
for k = 1:3
    figure(5+k); clf;
    do_surf(Xgrid, Tgrid, data_noisy{k}, ...
        sprintf('%s  (Cv = %.4f)', titles_noisy{k}, Cv_table(2,k)), cax, cmap);
    drawnow;
end

% -----------------------------------------------------------------------
% Figure 9 — eigenvalue spectra with analytical parabola overlay
% -----------------------------------------------------------------------
colors  = {[0.2 0.4 0.8], [1.0 0.5 0.0], [0.8 0.1 0.1]};
mkrs    = {'o', 's', 'd'};
methods = {'Exact DMD', 'fbDMD', 'TASC-DMD'};
omegas  = {{omega_edmd{1}, omega_fbdmd{1}, omega_tasc{1}}, ...
           {omega_edmd{2}, omega_fbdmd{2}, omega_tasc{2}}};
sp_ttls = {sprintf('No noise  (r = %d)', r1), ...
           sprintf('50%% noise  (r = %d)', r1)};

figure(9); clf;
tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

for di = 1:2
    nexttile;
    hold on; grid on;

    % Analytical parabola (ground truth)
    plot(Re_para, Im_range, 'k-', 'LineWidth', 2.5, ...
         'DisplayName', 'Analytical parabola (exact)');

    % Analytical eigenvalue points
    scatter(real(omega_exact), imag(omega_exact), 30, 'k', 'filled', ...
            'MarkerFaceAlpha', 0.25, 'DisplayName', 'Exact eigenvalues');

    % DMD eigenvalues with Cv in legend
    for mi = 1:3
        om = omegas{di}{mi};
        scatter(real(om), imag(om), 60, colors{mi}, ...
                'Marker', mkrs{mi}, 'LineWidth', 2, ...
                'DisplayName', sprintf('%s  (Cv=%.3f)', methods{mi}, Cv_table(di,mi)));
    end

    xlabel('Re($\omega_k$)', 'Interpreter','latex','FontSize',15,'FontWeight','bold');
    ylabel('Im($\omega_k$)', 'Interpreter','latex','FontSize',15,'FontWeight','bold');
    legend('Location','northeastoutside','FontSize', 9);
    title(sp_ttls{di}, 'FontSize', 12);
    xlim([-3  0.5]);
    ylim([-15  15]);
    axis square;
    set(gca, 'FontSize', 12, 'FontWeight', 'bold');
end

sgtitle({'Advection-Diffusion — Eigenvalue Spectrum vs Analytical Parabola', ...
         sprintf('Re(\\omega) = -(\\nu/U_0^2)\\cdotIm(\\omega)^2 = -%.4f\\cdotIm(\\omega)^2', nu/U0^2)}, ...
        'FontSize', 11, 'FontWeight', 'bold');
drawnow;

% -----------------------------------------------------------------------
% Figure 10 — Cv bar chart summary
% -----------------------------------------------------------------------
figure(10); clf;
b = bar(Cv_table, 0.75);
b(1).FaceColor = colors{1};
b(2).FaceColor = colors{2};
b(3).FaceColor = colors{3};
set(gca, 'XTickLabel', labels, 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Mean Cv  (global correlation)', 'FontSize', 13, 'FontWeight', 'bold');
title({'Advection-Diffusion — Reconstruction Quality', ...
       '(NOT in paper — additional benchmark)'}, ...
      'FontSize', 12, 'FontWeight', 'bold');
legend({'Exact DMD','fbDMD','TASC-DMD'}, 'Location', 'southwest', 'FontSize', 11);
ylim([0 1.05]);  grid on;
drawnow;
