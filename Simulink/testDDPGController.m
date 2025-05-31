function testDDPGController()
rng(30)
% TESTDDPGCONTROLLER Test the trained DDPG agent for monopile control
%
% This function loads a trained DDPG agent and tests its performance
% on the monopile control task, providing visualization of results.

%% Load Trained Agent
try
    load('trainedDDPGAgent.mat', 'agent');
    fprintf('Loaded trained DDPG agent successfully.\n');
catch
    error('Could not load trained agent. Please run trainDDPGController first.');
end

%% Test Parameters
testDuration = 60; % seconds
timeStep = 0.01; % seconds
nSteps = round(testDuration / timeStep);

% Initialize arrays for logging
time = (0:nSteps-1) * timeStep;
tip_displacement = zeros(nSteps, 1);
tip_velocity = zeros(nSteps, 1);

wave_elevation = zeros(nSteps, 1);
control_actions = zeros(nSteps, 2);
guy_wire_forces = zeros(nSteps, 1);
rewards = zeros(nSteps, 1);

% Remove LSTM-specific observation buffer - MLP uses single observations
% No need for sequenceLength or obsBuffer

fprintf('Starting DDPG agent testing for %d seconds...\n', testDuration);

%% Run Test Simulation
for step = 1:nSteps
    % Generate wave elevation (more realistic irregular wave)
    wave_elevation(step) = generateWaveElevation(time(step));
    
    % Get current observation
    if step == 1
        % Initialize with zero displacement
        [tip_disp, tip_vel, ~, ~, guy_forces] = monopileSimulinkFunction(0, 0);
    else
        tip_disp = tip_displacement(step-1);
        tip_vel = tip_velocity(step-1);

        guy_forces = guy_wire_forces(step-1);
    end
    
    % Prepare current observation for MLP (no sequence needed)
    currentObs = [tip_disp;tip_vel; wave_elevation(step)];
    
    % Get action from agent (MLP can act immediately)
    action = getAction(agent, currentObs);
    
    % Debug: Display action information
    if step <= 3  % Only show for first few steps
        if isnumeric(action)
            fprintf('Step %d: action class = %s, size = [%s], value = %s\n', ...
                step, class(action), num2str(size(action)), mat2str(action));
        else
            fprintf('Step %d: action class = %s, size = [%s], value = (non-numeric)\n', ...
                step, class(action), num2str(size(action)));
        end
    end
    
    % Convert action to numeric if it's a cell array
    if iscell(action)
        action = cell2mat(action);
    end
    
    % Ensure action is a column vector
    action = action(:);
    
    % Handle case where agent only returns 1 action instead of 2
    if length(action) == 1
        % Duplicate the single action for both guy wires
        action = [action; action];
        if step <= 3
            fprintf('Warning: Agent returned only 1 action, duplicating for both guy wires\n');
        end
    elseif length(action) ~= 2
        % Default to zero actions if unexpected size
        action = [0; 0];
        if step <= 3
            fprintf('Warning: Agent returned %d actions, using [0; 0]\n', length(action));
        end
    end
    
    % Apply control action to simulation
    [tip_disp, tip_vel, tip_acc, base_moment, guy_forces] = ...
        monopileSimulinkFunction(action(1), action(2));
    
    % Log results
    tip_displacement(step) = tip_disp;
    control_actions(step, :) = action;
    guy_wire_forces(step) = guy_forces;
    
    % Calculate reward
    rewards(step) = -abs(tip_disp) - 0.5*(abs(action(1)) + abs(action(2)));
    
    % Progress indicator
    if mod(step, 1000) == 0
        fprintf('Completed %d/%d steps (%.1f%%)\n', step, nSteps, 100*step/nSteps);
    end
end

fprintf('Testing completed.\n');

%% Performance Analysis
fprintf('\n=== DDPG Agent Performance Analysis ===\n');
fprintf('Maximum displacement: %.4f m\n', max(abs(tip_displacement)));
fprintf('RMS displacement: %.4f m\n', rms(tip_displacement));
fprintf('Average reward: %.2f\n', mean(rewards));
fprintf('Control effort (RMS): %.4f m\n', rms(sqrt(sum(control_actions.^2, 2))));

%% Visualization
createPlots(time, tip_displacement, wave_elevation, control_actions, guy_wire_forces, rewards,testDuration);

%% Save Results
results = struct();
results.time = time;
results.tip_displacement = tip_displacement;
results.wave_elevation = wave_elevation;
results.control_actions = control_actions;
results.guy_wire_forces = guy_wire_forces;
results.rewards = rewards;

save('ddpg_test_results.mat', 'results');
fprintf('\nResults saved to ddpg_test_results.mat\n');

end

function wave_elev = generateWaveElevation(t)
% Generate realistic irregular wave elevation using multiple components

% Wave parameters
H_s = 5.0; % Significant wave height [m]
T_p = 10.0; % Peak period [s]
omega_p = 2*pi/T_p; % Peak frequency [rad/s]

% Multiple frequency components for irregular waves
freqs = [0.5, 0.7, 1.0, 1.3, 1.5, 2.0] * omega_p;
amps = [0.3, 0.8, 1.0, 0.6, 0.4, 0.2] * H_s/2;
phases = [0, pi/4, pi/2, 3*pi/4, pi, 5*pi/4];

% Sum of sinusoidal components
%wave_elev = sum(amps .* sin(freqs * t + phases));
wave_elev = (H_s/2 * sin(omega_p * t + phases(1)));

end

function createPlots(time, displacement, wave_elevation, control_actions, guy_wire_forces, rewards,testDuration)
% Create comprehensive plots of DDPG agent performance

figure('Position', [100, 100, 1200, 800]);

% Subplot 1: Displacement and Wave Elevation
subplot(3, 2, 1);
yyaxis left
plot(time, displacement, 'b-', 'LineWidth', 1.5);
ylabel('Tip Displacement [m]');
yyaxis right
plot(time, wave_elevation, 'k--', 'LineWidth', 1);
ylabel('Wave Elevation [m]');
xlabel('Time [s]');
title('Tower Response and Wave Excitation');
grid on;
legend('Tip Displacement', 'Wave Elevation', 'Location', 'best');

% Subplot 2: Control Actions
subplot(3, 2, 2);
plot(time, control_actions(:,1), 'r-', 'LineWidth', 1.5);
hold on;
plot(time, control_actions(:,2), 'g-', 'LineWidth', 1.5);
ylabel('Control Displacement [m]');
xlabel('Time [s]');
title('Guy Wire Control Actions');
grid on;
legend('Cable 1 (\Delta L_1)', 'Cable 2 (\Delta L_2)', 'Location', 'best');

% Subplot 3: Guy Wire Forces
subplot(3, 2, 3);
plot(time, guy_wire_forces/1000, 'k-', 'LineWidth', 1.5);
ylabel('Guy Wire Forces [kN]');
xlabel('Time [s]');
title('Total Guy Wire Forces');
grid on;

% Subplot 4: Rewards
subplot(3, 2, 4);
plot(time, rewards, 'm-', 'LineWidth', 1);
ylabel('Reward');
xlabel('Time [s]');
title('DDPG Agent Rewards');
grid on;

% Subplot 5: Displacement Statistics
subplot(3, 2, 5);
histogram(displacement, 50, 'Normalization', 'probability');
xlabel('Tip Displacement [m]');
ylabel('Probability');
title('Displacement Distribution');
grid on;

% Subplot 6: Control Effort
subplot(3, 2, 6);
control_effort = sqrt(sum(control_actions.^2, 2));
plot(time, control_effort, 'Color', [0.8, 0.4, 0.0], 'LineWidth', 1.5);
ylabel('Control Effort [m]');
xlabel('Time [s]');
title('Total Control Effort');
grid on;

sgtitle('DDPG Agent Performance Analysis', 'FontSize', 16, 'FontWeight', 'bold');

% Save figure
savefig('ddpg_performance_analysis.fig');
print('ddpg_performance_analysis.png', '-dpng', '-r300');

end