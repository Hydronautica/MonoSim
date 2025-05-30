% MonopileSimulator - Main Script
% A modular implementation of the monopile simulator with clear separation of concerns
% This script serves as the "glue code" that calls various functions for:
% 1. Initialization: Read input parameters and create data structures
% 2. Simulation: Calculate external forces and solve Newmark integration
% 3. Visualization: Generate plots and animations

%% Clear workspace and set default parameters
clear; clc; close all;

% For reproducibility 
rng(30);
set(groot, 'DefaultAxesFontSize', 24); % Set default axes font size
set(groot, 'DefaultTextFontSize', 24); % Set default text font size

%% Initialize simulation parameters and data structures
params = initializeParameters();

%% Create mesh and assemble global matrices
addpath('FEM', 'Integration', 'soil', 'hydro', 'aero');
[K, M, C, elemZ, nNode, nDOF, free, fixed] = assembleMeshAndMatrices(params);

%% Initialize time integration variables
[U, V, Aacc, time, ramp] = initializeTimeIntegration(params, nDOF);

%% Compute effective stiffness matrix for Newmark integration
[Keff, a0c, a1c, a2c, a3c, a4c, a5c] = computeNewmarkConstants(params, K, M, C, free);

%% Initialize history arrays for forces and other quantities
[FHyd_history, FAero_history, Fsoil_history] = initializeForceHistory(params, nDOF);

%% Initialize soil state variables
[soilStates, thetaStates] = initializeSoilStates(params, nNode);

%% Compute initial acceleration
% Set step counter in params for force history tracking
params.step = 1;

% Compute initial forces
[F0, FHyd_history] = computeMorisonForces(U(:,1), V(:,1), time(1), elemZ, params, FHyd_history);
[Fm_air, FAero_history] = computeAerodynamicForces(time(1), elemZ, params, FAero_history);
F_const = computeConstantForces(params, nDOF, time(1));
[F_py, soilStates, thetaStates] = computeSoilForces(U(:,1), V(:,1), elemZ, params, soilStates, thetaStates);

% Add guy wire forces if enabled
F_guy = zeros(nDOF, 1);
if params.guywires
    F_guy = computeGuyWireForces(U(:,1), V(:,1), params);
    F_py = F_py + F_guy;
end

% Apply ramp to initial forces
ramp_i = ramp(1);
F0 = (F0 + F_const + Fm_air + F_py) * ramp_i;

% Compute initial acceleration
Aacc(free, 1) = M(free, free) \ (F0(free) - C(free, free) * V(free, 1) - K(free, free) * U(free, 1));

%% Initialize second-order wave elevation if needed
eta2 = 0;
eta2_history = zeros(size(time));
if params.secondOrder && params.irregular
    eta2 = secondOrderElevation(time(1), params);
    eta2_history(1) = eta2;
    params.eta2 = eta2;
end

%% Main simulation loop
disp('Starting main simulation loop...');
tic;

for i = 1:params.nSteps
    % Update step counter in params
    params.step = i + 1;
    
    % Print progress
    if mod(i, floor(params.nSteps/10)) == 0
        fprintf('Progress: %.1f%% (t = %.2f s)\n', 100*i/params.nSteps, time(i+1));
    end
    
    % Get current ramp factor
    ramp_i = ramp(i+1);
    
    % Calculate second-order wave elevation if needed
    if params.secondOrder && params.irregular
        eta2 = secondOrderElevation(time(i+1), params);
        eta2_history(i+1) = eta2;
        params.eta2 = eta2;
    end
    
    % Get external forces at current time step
    [F_hydro, FHyd_history] = computeMorisonForces(U(:,i), V(:,i), time(i+1), elemZ, params, FHyd_history);
    [Fm_air, FAero_history] = computeAerodynamicForces(time(i+1), elemZ, params, FAero_history);
    F_const = computeConstantForces(params, nDOF, time(i+1));
    
    % Compute soil reaction forces
    [F_py, soilStates, thetaStates] = computeSoilForces(U(:,i), V(:,i), elemZ, params, soilStates, thetaStates);
    
    % Compute guy wire forces if enabled
    F_guy = zeros(nDOF, 1);
    if params.guywires
        F_guy = computeGuyWireForces(U(:,i), V(:,i), params);
        F_py = F_py + F_guy;
    end
    
    % Store soil forces history
    Fsoil_history(:, i+1) = F_py;
    
    % Total external force vector
    Fi_t = (F_hydro + F_const + Fm_air + F_py) * ramp_i;
    
    % Solve Newmark integration step
    [U(:,i+1), V(:,i+1), Aacc(:,i+1)] = solveNewmarkStep(U(:,i), V(:,i), Aacc(:,i), Fi_t, K, M, C, free, fixed, a0c, a1c, a2c, a3c, a4c, a5c, Keff, params.dt, params.gammaN);
    
    % Store force history
    FHyd_history(:,i+1) = F_hydro;
    FAero_history(:,i+1) = Fm_air + F_const;
end

toc;
disp('Simulation completed!');

%% Post-processing and visualization
disp('Generating plots and visualizations...');

% Plot tip displacement and bending stress
plotTipDisplacementAndStress(time, U, elemZ, params);

% Plot nodal acceleration
plotNodalAcceleration(time, Aacc, elemZ, params);

% Calculate and display natural frequencies
[freq_hz, alpha_ray, beta_ray] = calculateNaturalFrequencies(K, M, free);
%
% Plot surface elevation
plotSurfaceElevation(time, eta2_history, params);

% Generate 2D animation if requested
if params.video
    generate2DAnimation(time, U, FHyd_history, FAero_history, elemZ, params);
end

% Generate 3D animation if requested
if params.video
    generate3DAnimation(time, U, FHyd_history, FAero_history, Fsoil_history, elemZ, params);
end

disp('All visualizations completed!');