function [sys, x0, str, ts] = monopileSFunction(t, x, u, flag)
% MONOPILESFUNCTION Level-1 S-Function wrapper for monopile simulator
% This S-Function provides a direct interface to the monopile simulator
% for use in Simulink.
%
% Inputs:
%   u(1) - deltaL_1 (control displacement for right guy wire)
%   u(2) - deltaL_2 (control displacement for left guy wire)
%
% Outputs:
%   y(1) - tip_displacement
%   y(2) - tip_velocity  
%   y(3) - tip_acceleration
%   y(4) - base_moment
%   y(5) - guy_wire_forces

switch flag
    case 0 % Initialization
        [sys, x0, str, ts] = mdlInitializeSizes;
        
    case 3 % Output calculation
        sys = mdlOutputs(t, x, u);
        
    case 9 % Terminate
        sys = mdlTerminate(t, x, u);
        
    otherwise
        sys = [];
end

function [sys, x0, str, ts] = mdlInitializeSizes
% Initialize the S-Function sizes

sizes = simsizes;
sizes.NumContStates = 0;     % No continuous states
sizes.NumDiscStates = 0;     % No discrete states
sizes.NumOutputs = 5;        % 5 outputs
sizes.NumInputs = 2;         % 2 inputs
sizes.DirFeedthrough = 1;    % Direct feedthrough
sizes.NumSampleTimes = 1;    % One sample time

sys = simsizes(sizes);
x0 = [];                     % No initial states
str = [];                    % No state names
ts = [-1 0];                 % Inherited sample time

function sys = mdlOutputs(t, x, u)
% Calculate outputs

try
    % Get input values
    deltaL_1 = u(1);
    deltaL_2 = u(2);
    
    % Call the monopile simulator function
    [tip_displacement, tip_velocity, tip_acceleration, base_moment, guy_wire_forces] = ...
        monopileSimulinkFunction(deltaL_1, deltaL_2);
    
    % Set outputs
    sys = [tip_displacement; tip_velocity; tip_acceleration; base_moment; guy_wire_forces];
    
catch ME
    % Handle errors gracefully
    warning('Monopile S-Function error: %s', ME.message);
    
    % Set default outputs in case of error
    sys = [0; 0; 0; 0; 0];
end

function sys = mdlTerminate(t, x, u)
% Terminate function - called at end of simulation
% Clean up any resources if needed
sys = [];