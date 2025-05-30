function [Keff, a0c, a1c, a2c, a3c, a4c, a5c] = computeNewmarkConstants(params, K, M, C, free)
% COMPUTENEWMARKCONSTANTS Computes constants needed for Newmark-beta integration
%
% Inputs:
%   params - Structure containing simulation parameters
%   K      - Global stiffness matrix
%   M      - Global mass matrix
%   C      - Global damping matrix
%   free   - Indices of free DOFs
%
% Outputs:
%   Keff   - Effective stiffness matrix for Newmark integration
%   a0c    - Newmark constant
%   a1c    - Newmark constant
%   a2c    - Newmark constant
%   a3c    - Newmark constant
%   a4c    - Newmark constant
%   a5c    - Newmark constant

% Extract Newmark parameters
betaN = params.betaN;
gammaN = params.gammaN;
dt = params.dt;

% Compute Newmark constants
a0c = 1/(betaN * dt^2);
a1c = gammaN / (betaN * dt);
a2c = 1/(betaN * dt);
a3c = 1/(2*betaN) - 1;
a4c = gammaN / betaN - 1;
a5c = dt * (gammaN/(2*betaN) - 1);

% Compute effective stiffness matrix
Keff = K(free, free) + a0c*M(free, free) + a1c*C(free, free);

end