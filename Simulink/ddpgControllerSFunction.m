function [sys,x0,str,ts] = ddpgControllerSFunction(t,x,u,flag)
% DDPGCONTROLLERSFUNCTION S-Function for DDPG controller in Simulink
%
% This S-Function implements a DDPG controller that takes observations
% (tip displacement, wave elevation, reset signal) and outputs control
% actions (deltaL_1, deltaL_2) for guy wire control.
%
% Inputs:
%   u(1) - tip displacement [m]
%   u(2) - wave elevation [m] 
%   u(3) - reset signal (0 or 1)
%
% Outputs:
%   y(1) - deltaL_1 (guy wire 1 control) [m]
%   y(2) - deltaL_2 (guy wire 2 control) [m]

switch flag
    case 0 % Initialization
        [sys,x0,str,ts] = mdlInitializeSizes;
        
    case 1 % Derivatives (not used for discrete system)
        sys = [];
        
    case 2 % Update (discrete states)
        sys = mdlUpdate(t,x,u);
        
    case 3 % Outputs
        sys = mdlOutputs(t,x,u);
        
    case 4 % GetTimeOfNextVarHit (not used)
        sys = [];
        
    case 9 % Terminate
        sys = [];
        
    otherwise
        error(['Unhandled flag = ',num2str(flag)]);
end

%% Initialize Sizes
function [sys,x0,str,ts] = mdlInitializeSizes

sizes = simsizes;
sizes.NumContStates  = 0;   % No continuous states
sizes.NumDiscStates  = 202; % 200 for observation buffer + 2 for control outputs
sizes.NumOutputs     = 2;   % deltaL_1, deltaL_2
sizes.NumInputs      = 3;   % tip_displacement, wave_elevation, reset
sizes.DirFeedthrough = 1;   % Direct feedthrough
sizes.NumSampleTimes = 1;   % Single sample time

sys = simsizes(sizes);

% Initial discrete states (observation buffer + control outputs)
x0 = zeros(202, 1);

% No state names
str = [];

% Sample time (0.01 seconds to match simulation)
ts = [0.01 0];

%% Update States
function sys = mdlUpdate(t,x,u)

% Extract inputs
tip_displacement = u(1);
wave_elevation = u(2);
reset_signal = u(3);

% Extract observation buffer (first 200 states: 2 obs x 100 time steps)
obsBuffer = reshape(x(1:200), [2, 100]);

% Reset buffer if requested
if reset_signal > 0.5
    obsBuffer = zeros(2, 100);
else
    % Update observation buffer (sliding window)
    obsBuffer(:, 1:end-1) = obsBuffer(:, 2:end);
    obsBuffer(:, end) = [tip_displacement; wave_elevation];
end

% Store updated buffer back to states
sys = [reshape(obsBuffer, [200, 1]); x(201:202)];

%% Generate Outputs
function sys = mdlOutputs(t,x,u)

persistent agent isAgentLoaded

% Initialize agent on first call
if isempty(isAgentLoaded)
    isAgentLoaded = false;
    try
        % Load trained agent
        loadedData = load('trainedDDPGAgent.mat');
        agent = loadedData.agent;
        isAgentLoaded = true;
        fprintf('DDPG agent loaded successfully in S-Function.\n');
    catch ME
        warning('Could not load DDPG agent: %s. Using zero control.', ME.message);
        agent = [];
        isAgentLoaded = false;
    end
end

% Extract observation buffer from states
obsBuffer = reshape(x(1:200), [2, 100]);

% Get action from agent
if isAgentLoaded && ~isempty(agent)
    try
        % Prepare observation sequence for LSTM (transpose for correct format)
        obsSequence = obsBuffer'; % [100 x 2] format expected by LSTM
        
        % Get action from agent
        action = getAction(agent, obsSequence);
        deltaL_1 = action(1);
        deltaL_2 = action(2);
        
        % Clamp outputs to reasonable range
        deltaL_1 = max(-0.5, min(0.5, deltaL_1));
        deltaL_2 = max(-0.5, min(0.5, deltaL_2));
        
    catch ME
        % Fallback if agent prediction fails
        warning('Agent prediction failed: %s. Using zero control.', ME.message);
        deltaL_1 = 0;
        deltaL_2 = 0;
    end
else
    % No agent available, use zero control
    deltaL_1 = 0;
    deltaL_2 = 0;
end

% Output control actions
sys = [deltaL_1; deltaL_2];