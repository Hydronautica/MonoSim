function F_const = computeConstantForces(params, nDOF, t)
% COMPUTECONSTANTFORCES Computes constant and oscillatory tip forces
%
% Inputs:
%   params - Structure containing simulation parameters
%   nDOF   - Number of degrees of freedom
%   t      - Current time (optional, for oscillatory component)
%
% Outputs:
%   F_const - Vector of constant forces

% Initialize force vector
F_const = zeros(nDOF, 1);

% Get tip DOF
topDOF = 2*(params.nNode-1) + 1;

% Apply constant tip force
F_const(topDOF) = params.F_hub;

% Add oscillatory component if time is provided
if nargin > 2
    F_const(topDOF) = F_const(topDOF) + params.F_hub_oscill * cos(params.F_hub_freq * t);
end

end