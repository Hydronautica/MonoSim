function [FHyd_history, FAero_history, Fsoil_history] = initializeForceHistory(params, nDOF)
% INITIALIZEFORCEHISTORY Initialize arrays to store force history
%
% Inputs:
%   params - Structure containing simulation parameters
%   nDOF   - Number of degrees of freedom
%
% Outputs:
%   FHyd_history  - History of hydrodynamic forces
%   FAero_history - History of aerodynamic forces
%   Fsoil_history - History of soil reaction forces

% Extract parameters
nSteps = params.nSteps;

% Initialize force history arrays
FHyd_history = zeros(nDOF, nSteps+1);
FAero_history = zeros(nDOF, nSteps+1);
Fsoil_history = zeros(nDOF, nSteps+1);

end