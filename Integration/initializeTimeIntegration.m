function [U, V, Aacc, time, ramp] = initializeTimeIntegration(params, nDOF)
% INITIALIZETIMEINTEGRATION Initialize displacement, velocity, and acceleration arrays
%
% Inputs:
%   params - Structure containing simulation parameters
%   nDOF   - Number of degrees of freedom
%
% Outputs:
%   U      - Displacement array (nDOF x nSteps+1)
%   V      - Velocity array (nDOF x nSteps+1)
%   Aacc   - Acceleration array (nDOF x nSteps+1)
%   time   - Time vector
%   ramp   - Ramp function vector

% Extract parameters
nSteps = params.nSteps;
dt = params.dt;
t_total = params.t_total;
ramp_duration = params.ramp_duration;

% Initialize arrays
U = zeros(nDOF, nSteps+1);
V = zeros(nDOF, nSteps+1);
Aacc = zeros(nDOF, nSteps+1);

% Create time vector
time = (0:nSteps) * dt;

% Create ramp function
ramp = min(time / ramp_duration, 1);  % Linear ramp from 0 to 1 over ramp_duration seconds

end