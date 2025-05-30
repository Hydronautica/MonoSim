function [Fs, state] = pySimpleHysteretic(y, ydot, state, params, z)
% PYSIMPLEHYSTERETIC Computes soil reaction force using p-y hysteretic model
%
% Inputs:
%   y      - Current displacement
%   ydot   - Current velocity
%   state  - Current state of the p-y model
%   params - Structure containing simulation parameters
%   z      - Depth below mudline (positive downward)
%
% Outputs:
%   Fs     - Soil reaction force
%   state  - Updated state of the p-y model

% Extract parameters
linearPY = params.linearPY;
k_soil = params.k_soil;
Le = params.Le;
D_outer1 = params.D_outer1;

% If using linear p-y curves, compute simple spring force
if linearPY
    % Linear spring force
    Fs_spring = -k_soil * z * y * Le;
    
    % Add damping component (5% of critical damping)
    c_soil = 0.05 * 2 * sqrt(k_soil * z * Le);
    Fs_damping = -c_soil * ydot;
    
    % Total force
    Fs = Fs_spring + Fs_damping;
    return;
end

% Otherwise, use nonlinear hysteretic p-y model

% Extract state variables
y_prev = state.y_prev;
y_max = state.y_max;
y_min = state.y_min;
F_prev = state.F_prev;
loading = state.loading;

% Calculate ultimate soil resistance (API method)
z_abs = abs(z);  % Ensure positive depth
gamma_soil = 10000;  % Soil unit weight (N/m³)
A = min(0.9, 3 - 0.8*(z_abs/D_outer1));
pu = min(9*z_abs*gamma_soil*D_outer1, 3*z_abs*gamma_soil*D_outer1 + gamma_soil*z_abs^2*D_outer1*A);

% Initial stiffness
k_ini = k_soil * z_abs;

% Reference displacement at 50% of ultimate resistance
y50 = 0.01 * D_outer1;

% Check loading/unloading condition
if y > y_prev
    loading = true;
else
    loading = false;
end

% Update displacement extremes
if y > y_max
    y_max = y;
elseif y < y_min
    y_min = y;
end

% Calculate backbone curve force
F_backbone = sign(y) * pu * tanh(abs(k_ini * y / pu));

% Calculate force based on loading/unloading condition
if loading
    % Loading from min to max
    if y_min < 0
        % Coming from negative displacement
        F_range = F_backbone - (-pu * tanh(abs(k_ini * y_min / pu)));
        y_range = y - y_min;
        F_s = F_prev + F_range * (y - y_prev) / y_range;
    else
        % Loading from zero or positive
        F_s = F_backbone;
    end
else
    % Unloading from max to min
    if y_max > 0
        % Coming from positive displacement
        F_range = F_backbone - (pu * tanh(abs(k_ini * y_max / pu)));
        y_range = y - y_max;
        F_s = F_prev + F_range * (y - y_prev) / y_range;
    else
        % Unloading from zero or negative
        F_s = F_backbone;
    end
end

% Add damping component (5% of critical damping)
c_soil = 0.05 * 2 * sqrt(k_soil * z_abs);
F_damping = -c_soil * ydot;

% Total force
Fs = F_s * Le + F_damping;

% Update state
state.y_prev = y;
state.y_max = y_max;
state.y_min = y_min;
state.F_prev = F_s;
state.loading = loading;

end