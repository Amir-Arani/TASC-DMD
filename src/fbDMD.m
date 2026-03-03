function [nt2, Phi_fb, Lambda_fb, omega_fb, time_dynamics_fb, X_fbdmd, Cv] = fbDMD(X, dt, p, tde, r)
% FBDMD  Forward-Backward Dynamic Mode Decomposition with optional
%        time-delay embedding.
%
%   [nt2, Phi_fb, Lambda_fb, omega_fb, time_dynamics_fb, X_fbdmd, Cv] = ...
%       fbDMD(X, dt, p, tde, r)
%
%   Computes fbDMD by combining the forward and backward linear operators
%   via a geometric mean (matrix square root), yielding eigenvalues that
%   are symmetric about the imaginary axis — improving noise robustness.
%
%   Inputs:
%       X    - Data matrix  [n x m]  (n: spatial DOFs, m: snapshots)
%       dt   - Time step between snapshots
%       p    - Temporal subsampling factor (effective dt becomes p*dt)
%       tde  - Number of time delays (1, 2, or 3)
%       r    - Rank truncation for the initial SVD
%
%   Outputs:
%       nt2              - Number of time points after embedding
%       Phi_fb           - fbDMD modes  [n*tde x r]
%       Lambda_fb        - Discrete-time eigenvalues  [r x 1]
%       omega_fb         - Continuous-time eigenvalues  [r x 1]
%       time_dynamics_fb - Temporal coefficients  [r x nt2]
%       X_fbdmd          - fbDMD reconstruction  [n*tde x nt2]
%       Cv               - Global correlation coefficient per snapshot  [nt2 x 1]
%
%   Reference:
%       [TASC DMD Paper — citation to be added after publication]

%       Dawson et al. (2016), "Characterizing and correcting for the effect
%       of sensor noise in the dynamic mode decomposition."
%
%   Last updated: 10/21/2024 by AHGA

% -------------------------------------------------------------------------
% Build time-delay embedded data matrix
% -------------------------------------------------------------------------
x0 = X(:, 1:p:end);
nt1 = size(x0, 2);

x = [];
for i = 1 : nt1 - tde + 1
    switch tde
        case 1,  xx = x0(:, i);
        case 2,  xx = [x0(:, i); x0(:, i+1)];
        case 3,  xx = [x0(:, i); x0(:, i+1); x0(:, i+2)];
    end
    x = [x, xx]; %#ok<AGROW>
end

nt2 = size(x, 2);
X1 = x(:, 1:end-1);
X2 = x(:, 2:end);

% -------------------------------------------------------------------------
% Initial rank-r truncation
% -------------------------------------------------------------------------
[U_r, S_r, V_r] = svd(X1, 'econ');
Ur = U_r(:, 1:r);
Sr = S_r(1:r, 1:r);
Vr = V_r(:, 1:r);

XX = Ur' * X1;   % Projected forward snapshots
YY = Ur' * X2;   % Projected backward snapshots

% -------------------------------------------------------------------------
% Forward operator (KF) and backward operator (KB) in the projected space
% -------------------------------------------------------------------------
[Uf, Sf, Vf] = svd(XX, 'econ');
[Ub, Sb, Vb] = svd(YY, 'econ');

KFtilde = Uf' * YY * Vf / Sf;
KBtilde = Ub' * XX * Vb / Sb;

SF = YY * Vf / Sf;
SB = XX * Vb / Sb;

KF = SF * KFtilde * pinv(SF);
KB = SB * KBtilde * pinv(SB);

% Geometric mean: A_fb = sqrt(KF * KB^{-1})
Atilde_fb = sqrtm(KF / KB);
[W_fb, D_fb] = eig(Atilde_fb);

% -------------------------------------------------------------------------
% fbDMD modes, eigenvalues, and temporal reconstruction
% -------------------------------------------------------------------------
Phi_fb   = X2 * (Vr / Sr) * W_fb;
Lambda_fb = diag(D_fb);
omega_fb  = log(Lambda_fb) / (p * dt);

b_fb = Phi_fb \ X1(:, 1);

time_dynamics_fb = zeros(r, nt2);
for i = 1 : nt2
    time_dynamics_fb(:, i) = b_fb .* exp(omega_fb * (i-1) * (p*dt));
end
X_fbdmd = Phi_fb * time_dynamics_fb;

% -------------------------------------------------------------------------
% Global correlation coefficient Cv
% -------------------------------------------------------------------------
Cv_cmplx = zeros(nt2, 1);
Cv       = zeros(nt2, 1);
for j = 1 : nt2
    Xn_org = normalize(x(1:end/tde, j), 'norm');
    Xn_dmd = normalize(X_fbdmd(1:end/tde, j), 'norm');
    Cv_cmplx(j) = dot(Xn_org, Xn_dmd);
    Cv(j)       = abs(Cv_cmplx(j));
end

end
