function [nt2, Phi, lambda, omega, time_dynamics, X_dmd, Cv] = DMD_exact(X, dt, p, tde, r, var)
% DMD_EXACT  Exact (or standard) Dynamic Mode Decomposition with optional
%            time-delay embedding.
%
%   [nt2, Phi, lambda, omega, time_dynamics, X_dmd, Cv] = ...
%       DMD_exact(X, dt, p, tde, r, var)
%
%   Inputs:
%       X    - Data matrix  [n x m]  (n: spatial DOFs, m: snapshots)
%       dt   - Time step between snapshots
%       p    - Temporal subsampling factor (effective dt becomes p*dt)
%       tde  - Number of time delays (1, 2, or 3)
%       r    - Number of DMD modes (rank truncation)
%       var  - DMD variant: 0 = Exact DMD,  1 = Standard DMD
%
%   Outputs:
%       nt2            - Number of time points after embedding
%       Phi            - DMD modes  [n*tde x r]
%       lambda         - Discrete-time eigenvalues  [r x 1]
%       omega          - Continuous-time eigenvalues  [r x 1]
%       time_dynamics  - Temporal coefficients  [r x nt2]
%       X_dmd          - DMD reconstruction  [n*tde x nt2]
%       Cv             - Global correlation coefficient per snapshot  [nt2 x 1]
%
%   Reference:
%       [TASC DMD Paper — citation to be added after publication]

%
%   Last updated: 10/21/2024 by AHGA

% -------------------------------------------------------------------------
% Build time-delay embedded data matrix
% -------------------------------------------------------------------------
x0 = X(:, 1:p:end);          % Subsample at interval p
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
% SVD and reduced-order linear operator
% -------------------------------------------------------------------------
[U_r, S_r, V_r] = svd(X1, 'econ');
Ur = U_r(:, 1:r);
Sr = S_r(1:r, 1:r);
Vr = V_r(:, 1:r);

Atilde = Ur' * X2 * Vr / Sr;
[W, D] = eig(Atilde);

% -------------------------------------------------------------------------
% DMD modes
% -------------------------------------------------------------------------
switch var
    case 0
        Phi = X2 * (Vr / Sr) * W;   % Exact DMD modes
    case 1
        Phi = Ur * W;                % Standard DMD modes
end

lambda = diag(D);
omega  = log(lambda) / (p * dt);

% -------------------------------------------------------------------------
% Amplitude coefficients (least-squares fit to initial condition)
% -------------------------------------------------------------------------
b = Phi \ X1(:, 1);

% -------------------------------------------------------------------------
% Reconstruct temporal dynamics and full field
% -------------------------------------------------------------------------
time_dynamics = zeros(r, nt2);
for v = 1 : nt2
    time_dynamics(:, v) = b .* exp(omega * (v-1) * (p*dt));
end
X_dmd = Phi * time_dynamics;

% -------------------------------------------------------------------------
% Global correlation coefficient Cv
% -------------------------------------------------------------------------
Cv_cmplx = zeros(nt2, 1);
Cv       = zeros(nt2, 1);
for j = 1 : nt2
    Xn_org = normalize(x(1:end/tde, j), 'norm');
    Xn_dmd = normalize(X_dmd(1:end/tde, j), 'norm');
    Cv_cmplx(j) = dot(Xn_org, Xn_dmd);
    Cv(j)       = abs(Cv_cmplx(j));
end

end
