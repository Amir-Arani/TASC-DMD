function rhs = nls_rhs(t, psit, dummy, k)
% NLS_RHS  Right-hand side of the Nonlinear Schrödinger (NLS) equation
%          in Fourier space, for use with ODE solvers (e.g., ode45).
%
%   Equation:  i*psi_t + (1/2)*psi_xx + |psi|^2 * psi = 0
%
%   Usage:
%       rhs = nls_rhs(t, psit, [], k)
%
%   Inputs:
%       t     - Current time (unused, required by ODE solver interface)
%       psit  - Fourier coefficients of psi at current time  [n x 1]
%       dummy - Unused placeholder (required by ode45 interface)
%       k     - Wave number vector  [n x 1]
%
%   Output:
%       rhs   - Time derivative of psit in Fourier space  [n x 1]
%
%   Note:
%       The split-step Fourier method is used:
%         - Linear term:    -(i/2) * k^2 * psit
%         - Nonlinear term: +i * FFT(|psi|^2 * psi)

psi = ifft(psit);
rhs = -(1i/2) .* (k.^2) .* psit + 1i .* fft((abs(psi).^2) .* psi);

end
