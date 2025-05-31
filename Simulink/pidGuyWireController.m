function [control_guy1,time,displacement, time2, displacementBaseline] = pidGuyWireController()
% PIDGUYWIRECONTROLLER PID controller for individual guy wire control
%
% This function implements PID controllers for each guy wire to minimize
% tower top displacement. Each guy wire has its own PID controller that
% responds to the tower displacement and velocity.

%% Clear workspace and set random seed
clear; clc; close all;
rng(44); % For reproducibility
clear monopileSimulinkFunction;

%% PID Controller Parameters
% Guy Wire 1 PID parameters
Kp1 = 0.3;    % Proportional gain
Ki1 = 0.05;    % Integral gain  
Kd1 = 0.25;   % Derivative gain

% Guy Wire 2 PID parameters
Kp2 = 0.3;    % Proportional gain
Ki2 = 0.05;    % Integral gain
Kd2 = 0.25;   % Derivative gain

%% Simulation Parameters
dt = 0.005;           % Time step [s]
simTime = 3600;       % Total simulation time [s]
nSteps = simTime/dt; % Number of simulation steps

% Target displacement (ideally zero)
target_displacement = 0;

%% Initialize variables
time = (0:nSteps-1) * dt;
time2 = 0:0.01:199.99;

displacement = zeros(1, nSteps);
velocity = zeros(1, nSteps);
control_guy1 = zeros(1, nSteps);
control_guy2 = zeros(1, nSteps);

% PID controller states
error_guy1 = zeros(1, nSteps);
error_guy2 = zeros(1, nSteps);
integral_guy1 = 0;
integral_guy2 = 0;
prev_error_guy1 = 0;
prev_error_guy2 = 0;

% Wave parameters for disturbance
wave_amplitude = 2.5;  % Wave amplitude [m]
wave_period = 10;      % Wave period [s]
wave_frequency = 2*pi/wave_period;

%% Main simulation loop
fprintf('Starting PID guy wire control simulation...\n');
fprintf('Simulation time: %.1f seconds\n', simTime);
fprintf('Time step: %.3f seconds\n', dt);

for i = 1:nSteps
    current_time = time(i);
    
    % Generate wave disturbance
    wave_elevation = wave_amplitude * sin(wave_frequency * current_time);
    
    % Calculate error for each guy wire (both respond to same displacement)
    error_guy1(i) = -(target_displacement - displacement(i));
    error_guy2(i) =(target_displacement - displacement(i));
    
    % PID Controller for Guy Wire 1
    integral_guy1 = integral_guy1 + error_guy1(i) * dt;
    derivative_guy1 = (error_guy1(i) - prev_error_guy1) / dt;
    control_guy1(i) = Kp1 * error_guy1(i) + Ki1 * integral_guy1 + Kd1 * derivative_guy1;
    
    % PID Controller for Guy Wire 2
    integral_guy2 = integral_guy2 + error_guy2(i) * dt;
    derivative_guy2 = (error_guy2(i) - prev_error_guy2) / dt;
    control_guy2(i) = Kp2 * error_guy2(i) + Ki2 * integral_guy2 + Kd2 * derivative_guy2;
    
    % Apply control limits (guy wire displacement limits)
    control_guy1(i) = max(-2, min(2, control_guy1(i)));
    control_guy2(i) =  max(-2, min(2, control_guy2(i)));
    
    % Call monopile simulation function
    try
        [tip_displacement, tip_velocity, ~, ~, guy_wire_forces] = ...
            monopileSimulinkFunction(1*control_guy1(i),1*control_guy2(i));
        
        % Store results
        if i < nSteps
            displacement(i+1) = tip_displacement;
            velocity(i+1) = tip_velocity;
        end
        
    catch ME
        fprintf('Error at step %d: %s\n', i, ME.message);
        break;
    end
    
    % Update previous errors
    prev_error_guy1 = error_guy1(i);
    prev_error_guy2 = error_guy2(i);
    
    % Progress indicator
    if mod(i, round(nSteps/10)) == 0
        fprintf('Progress: %.0f%%\n', (i/nSteps)*100);
    end
end

%% Performance Analysis
fprintf('\nSimulation completed!\n');

% Calculate performance metrics
rms_displacement = sqrt(mean(displacement.^2));
max_displacement = max(abs(displacement));
rms_control_guy1 = sqrt(mean(control_guy1.^2));
rms_control_guy2 = sqrt(mean(control_guy2.^2));
max_control_guy1 = max(abs(control_guy1));
max_control_guy2 = max(abs(control_guy2));

fprintf('\nPerformance Metrics:\n');
fprintf('RMS Displacement: %.4f m\n', rms_displacement);
fprintf('Max Displacement: %.4f m\n', max_displacement);
fprintf('RMS Control Guy Wire 1: %.4f m\n', rms_control_guy1);
fprintf('RMS Control Guy Wire 2: %.4f m\n', rms_control_guy2);
fprintf('Max Control Guy Wire 1: %.4f m\n', max_control_guy1);
fprintf('Max Control Guy Wire 2: %.4f m\n', max_control_guy2);

%% Plotting Results
figure('Position', [100, 100, 1200, 800]);
load displacementBaseline.mat
load velocityBaseline.mat
% Displacement plot
subplot(1, 2, 1);
plot(time, displacement, 'b-', 'LineWidth', 3.5);
hold on
plot(time, displacementBaseline, 'LineWidth', 3.5);
axis([0 50 -1.6 1.5])
grid on;
xlabel('Time [s]');
ylabel('Displacement [m]');
title('Tower Top Displacement');
legend 'Control' 'No Control'
% % Velocity plot
% subplot(3, 2, 2);
% plot(time, velocity, 'r-', 'LineWidth', 1.5);
% hold on
% plot(time, velocityBaseline, 'r-', 'LineWidth', 1.5);
% 
% grid on;
% xlabel('Time [s]');
% ylabel('Velocity [m/s]');
% title('Tower Top Velocity');
% legend 'Control' 'No Control'

% Control signals
subplot(1, 2, 2);
plot(time, control_guy1, '-', 'LineWidth', 3.5);
hold on;
plot(time, control_guy2, '-', 'LineWidth', 3.5);
grid on;
xlabel('Time [s]');
ylabel('Control Signal [m]');
title('Guy Wire Control Signals');
legend('Guy Wire 1', 'Guy Wire 2', 'Location', 'best');
axis([0 50 -0.3 0.3])

% % Error signals
% subplot(3, 2, 4);
% plot(time, error_guy1, 'c-', 'LineWidth', 1.5);
% hold on;
% plot(time, error_guy2, 'k-', 'LineWidth', 1.5);
% grid on;
% xlabel('Time [s]');
% ylabel('Error [m]');
% title('PID Controller Errors');
% legend('Guy Wire 1 Error', 'Guy Wire 2 Error', 'Location', 'best');
% 
% % Control effort comparison
% subplot(3, 2, 5);
% control_effort = abs(control_guy1) + abs(control_guy2);
% plot(time, control_effort, 'r-', 'LineWidth', 1.5);
% grid on;
% xlabel('Time [s]');
% ylabel('Total Control Effort [m]');
% title('Total Control Effort');
% 
% % Phase plot (displacement vs velocity)
% subplot(3, 2, 6);
% plot(displacement, velocity, 'b-', 'LineWidth', 1);
% grid on;
% xlabel('Displacement [m]');
% ylabel('Velocity [m/s]');
% title('Phase Plot (Displacement vs Velocity)');
% 
% sgtitle('PID Guy Wire Controller Performance Analysis');

%% Save results
save('pidControllerResults.mat', 'time', 'displacement', 'velocity', ...
     'control_guy1', 'control_guy2', 'error_guy1', 'error_guy2', ...
     'Kp1', 'Ki1', 'Kd1', 'Kp2', 'Ki2', 'Kd2');

fprintf('\nResults saved to pidControllerResults.mat\n');
fprintf('PID controller analysis complete!\n');

end

%% Helper function for PID tuning
function tunePIDParameters()
% TUNEPIDPARAMETERS Interactive PID parameter tuning
%
% This function provides a framework for tuning PID parameters
% You can modify the parameters and re-run the simulation

fprintf('\nPID Parameter Tuning Guidelines:\n');
fprintf('Kp (Proportional): Controls response speed and steady-state error\n');
fprintf('  - Increase Kp to reduce steady-state error\n');
fprintf('  - Too high Kp causes oscillations\n');
fprintf('\n');
fprintf('Ki (Integral): Eliminates steady-state error\n');
fprintf('  - Increase Ki to eliminate steady-state error\n');
fprintf('  - Too high Ki causes instability\n');
fprintf('\n');
fprintf('Kd (Derivative): Reduces overshoot and improves stability\n');
fprintf('  - Increase Kd to reduce overshoot\n');
fprintf('  - Too high Kd amplifies noise\n');
fprintf('\n');
fprintf('Recommended starting values:\n');
fprintf('Kp = 0.1 to 1.0\n');
fprintf('Ki = 0.01 to 0.5\n');
fprintf('Kd = 0.001 to 0.1\n');

end