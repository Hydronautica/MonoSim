function [tip_displacement, tip_velocity, tip_acceleration, base_moment, guy_wire_forces] = monopileSimulinkFunction(deltaL_1, deltaL_2)
% MONOPILESIMULINKFUNCTION Full dynamics monopile simulation for Simulink
% 
% This function implements the complete monopile simulator with full time-domain
% dynamics, optimized for Simulink with persistent variables to avoid reinitialization.
%
% Inputs:
%   deltaL_1 - Control displacement for guy wire 1 (right) [m]
%   deltaL_2 - Control displacement for guy wire 2 (left) [m]
%
% Outputs:
%   tip_displacement - Lateral displacement at tower tip [m]
%   tip_velocity     - Lateral velocity at tower tip [m/s]
%   tip_acceleration - Lateral acceleration at tower tip [m/s^2]
%   base_moment      - Bending moment at mudline [N*m]
%   guy_wire_forces  - Forces from both guy wires [N]

% Use persistent variables to avoid reinitialization
persistent params K M C elemZ nNode nDOF free fixed initialized
persistent U V Aacc time ramp Keff a0c a1c a2c a3c a4c a5c
persistent FHyd_history FAero_history Fsoil_history soilStates thetaStates
persistent waveVel waveAcc eta2_history current_step

% Initialize only once
if isempty(initialized)
    % Add paths for required functions
    addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'FEM'));
    addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'Integration'));
    addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'soil'));
    addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'hydro'));
    addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'aero'));
    
    % Initialize simulation parameters
    params = initializeParameters();
    params.guywires = true;
    
    % Create mesh and assemble global matrices
    [K, M, C, elemZ, nNode, nDOF, free, fixed] = assembleMeshAndMatrices(params);
    
    % Initialize time integration variables
    [U, V, Aacc, time, ramp] = initializeTimeIntegration(params, nDOF);
    
    % Compute effective stiffness matrix for Newmark integration
    [Keff, a0c, a1c, a2c, a3c, a4c, a5c] = computeNewmarkConstants(params, K, M, C, free);
    
    % Initialize history arrays for forces
    [FHyd_history, FAero_history, Fsoil_history] = initializeForceHistory(params, nDOF);
    
    % Initialize soil state variables
    [soilStates, thetaStates] = initializeSoilStates(params, nNode);
    
    % Pre-compute wave kinematics for all time steps
    [waveVel, waveAcc, eta2_history] = precomputeWaveKinematics(time, elemZ, params);
    
    % Initialize step counter
    current_step = 1;
    
    % Compute initial acceleration
    params.step = 1;
    
    % Compute initial forces
    [F0, FHyd_history] = computeMorisonForces(waveVel(:,1), waveAcc(:,1), time(1), elemZ, params, FHyd_history, waveVel(:,1), waveAcc(:,1));
    [Fm_air, FAero_history] = computeAerodynamicForces(time(1), elemZ, params, FAero_history);
    F_const = computeConstantForces(params, nDOF, time(1));
    [F_py, soilStates, thetaStates] = computeSoilForces(U(:,1), V(:,1), elemZ, params, soilStates, thetaStates);
    
    % Add guy wire forces
    F_guy = zeros(nDOF, 1);
    if params.guywires
        F_guy = computeGuyWireForces(U(:,1), V(:,1), params);
    end
    
    % Apply ramp to initial forces
    ramp_i = ramp(1);
    F0 =  (F0 + F_const + Fm_air + F_py + F_guy) * ramp_i;
    
    % Compute initial acceleration
    Aacc(free, 1) = M(free, free) \ (F0(free) - C(free, free) * V(free, 1) - K(free, free) * U(free, 1));
    
    % Initialize second-order wave elevation if needed
    if params.secondOrder && params.irregular
        eta2 = secondOrderElevation(time(1), params);
        eta2_history(1) = eta2;
        params.eta2 = eta2;
    end
    
    initialized = true;
end

% Update guy wire parameters based on control inputs
params.L0_guyR = params.L0_guyR - deltaL_1;
params.L0_guyL = params.L0_guyL - deltaL_2;

% Advance simulation by one time step
current_step = current_step + 1;
i = current_step - 1; % Array index for previous step

% Handle cycling through pre-computed data when reaching the end
wave_step = mod(current_step - 1, params.nSteps) + 1;
wave_step_prev = mod(i - 1, params.nSteps) + 1;

% Update step counter
params.step = current_step;

% Get current ramp factor (use final ramp value after ramp period)
if wave_step <= length(ramp)
    ramp_i = ramp(wave_step);
else
    ramp_i = 1.0; % Full force after ramp period
end

% Update eta2 if needed (cycle through pre-computed values)
if params.secondOrder && params.irregular
    params.eta2 = eta2_history(wave_step);
end

% Expand arrays if needed to accommodate longer simulation
if current_step > size(U, 2)
    % Expand state arrays by doubling size
    new_size = size(U, 2) * 2;
    U(:, end+1:new_size) = 0;
    V(:, end+1:new_size) = 0;
    Aacc(:, end+1:new_size) = 0;
    FHyd_history(:, end+1:new_size) = 0;
    FAero_history(:, end+1:new_size) = 0;
    Fsoil_history(:, end+1:new_size) = 0;
end

% Get external forces at current time step (using cyclic wave data)
[F_hydro, FHyd_history] = computeMorisonForces(waveVel(:,wave_step_prev), waveAcc(:,wave_step_prev), time(wave_step), elemZ, params, FHyd_history, waveVel(:,wave_step), waveAcc(:,wave_step));
[Fm_air, FAero_history] = computeAerodynamicForces(time(wave_step), elemZ, params, FAero_history);
F_const = computeConstantForces(params, nDOF, time(wave_step));

% Compute soil reaction forces
[F_py, soilStates, thetaStates] = computeSoilForces(U(:,i), V(:,i), elemZ, params, soilStates, thetaStates);

% Compute guy wire forces
F_guy = zeros(nDOF, 1);
if params.guywires
    F_guy = computeGuyWireForces(U(:,i), V(:,i), params);
    F_py = F_py + F_guy;
end

% Store soil forces history
Fsoil_history(:, current_step) = F_py;

% Total external force vector
Fi_t = (F_hydro + F_const + Fm_air + F_py) * ramp_i;

% Solve Newmark integration step
[U(:,current_step), V(:,current_step), Aacc(:,current_step)] = solveNewmarkStep(U(:,i), V(:,i), Aacc(:,i), Fi_t, K, M, C, free, fixed, a0c, a1c, a2c, a3c, a4c, a5c, Keff, params.dt, params.gammaN);

% Store force history
FHyd_history(:,current_step) = F_hydro;
FAero_history(:,current_step) = Fm_air + F_const;

%% Extract outputs
% Find tip node (highest elevation)
[~, tip_node_idx] = max(elemZ);
tip_dof = 2*(tip_node_idx-1) + 1; % Horizontal DOF for tip node

% Get current values
tip_displacement = U(tip_dof, current_step);
tip_velocity = V(tip_dof, current_step);
tip_acceleration = Aacc(tip_dof, current_step);

% Calculate base moment (bending moment at mudline)
% Find node closest to mudline (z = 0)
[~, mudline_node_idx] = min(abs(elemZ));
mudline_dof_rot = 2*mudline_node_idx; % Rotational DOF at mudline
base_moment = K(mudline_dof_rot, :) * U(:, current_step); % Moment from stiffness matrix

% Calculate total guy wire forces
if params.guywires
    F_guy_final = computeGuyWireForces(U(:,current_step), V(:,current_step), params);
    guy_wire_forces = norm(F_guy_final); % Magnitude of guy wire force vector
else
    guy_wire_forces = 0;
end

end