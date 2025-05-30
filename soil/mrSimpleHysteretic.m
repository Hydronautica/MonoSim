function [Ms, state] = mrSimpleHysteretic(theta, thetadot, state, params, z)
% MRSIMPLEHYSTERETIC Computes soil moment reaction using m-theta hysteretic model
%
% Inputs:
%   theta     - Current rotation
%   thetadot  - Current rotational velocity
%   state     - Current state of the m-theta model
%   params    - Structure containing simulation parameters
%   z         - Depth below mudline (positive downward)
%
% Outputs:
%   Ms        - Soil reaction moment
%   state     - Updated state of the m-theta model

% Extract parameters
linearPY = params.linearPY;
k_soil_rot = params.k_soil_rot;
Le = params.Le;
D_outer1 = params.D_outer1;

% If using linear m-theta curves, compute simple spring moment
if linearPY
    % Linear spring moment
    Ms_spring = -k_soil_rot * z^2 * theta * Le;
    
    % Add damping component (5% of critical damping)
    c_soil_rot = 0.05 * 2 * sqrt(k_soil_rot * z^2 * Le);
    Ms_damping = -c_soil_rot * thetadot;
    
    % Total moment
    Ms = Ms_spring + Ms_damping;
    return;
end

% Otherwise, use nonlinear hysteretic m-theta model

% Extract state variables
theta_prev = state.theta_prev;
theta_max = state.theta_max;
theta_min = state.theta_min;
M_prev = state.M_prev;
loading = state.loading;

% Calculate ultimate soil moment resistance
z_abs = abs(z);  % Ensure positive depth
gamma_soil = 10000;  % Soil unit weight (N/m³)
mu = 0.4 * gamma_soil * z_abs^3 * D_outer1;  % Ultimate moment capacity

% Initial rotational stiffness
k_ini_rot = k_soil_rot * z_abs^2;

% Reference rotation at 50% of ultimate resistance
theta50 = 0.01;  % radians

% Check loading/unloading condition
if theta > theta_prev
    loading = true;
else
    loading = false;
end

% Update rotation extremes
if theta > theta_max
    theta_max = theta;
elseif theta < theta_min
    theta_min = theta;
end

% Calculate backbone curve moment
M_backbone = sign(theta) * mu * tanh(abs(k_ini_rot * theta / mu));

% Calculate moment based on loading/unloading condition
if loading
    % Loading from min to max
    if theta_min < 0
        % Coming from negative rotation
        M_range = M_backbone - (-mu * tanh(abs(k_ini_rot * theta_min / mu)));
        theta_range = theta - theta_min;
        M_s = M_prev + M_range * (theta - theta_prev) / theta_range;
    else
        % Loading from zero or positive
        M_s = M_backbone;
    end
else
    % Unloading from max to min
    if theta_max > 0
        % Coming from positive rotation
        M_range = M_backbone - (mu * tanh(abs(k_ini_rot * theta_max / mu)));
        theta_range = theta - theta_max;
        M_s = M_prev + M_range * (theta - theta_prev) / theta_range;
    else
        % Unloading from zero or negative
        M_s = M_backbone;
    end
end

% Add damping component (5% of critical damping)
c_soil_rot = 0.05 * 2 * sqrt(k_soil_rot * z_abs^2);
M_damping = -c_soil_rot * thetadot;

% Total moment
Ms = M_s * Le + M_damping;

% Update state
state.theta_prev = theta;
state.theta_max = theta_max;
state.theta_min = theta_min;
state.M_prev = M_s;
state.loading = loading;

end