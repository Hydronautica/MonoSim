function trainDDPGController(agent)
% Set random seed for reproducibility
seed = 31;
rng(seed, 'twister');  % Set MATLAB's random number generator

% For deep learning reproducibility (if using Deep Learning Toolbox)
if exist('rng', 'builtin')
    rng(seed, 'twister');
end

% Set environment variable for additional reproducibility
setenv('PYTHONSEED', num2str(seed));
% TRAINDDPGCONTROLLER Train DDPG agent with LSTM networks for monopile control
%
% This function implements a Deep Deterministic Policy Gradient (DDPG) agent
% with LSTM networks to learn optimal guy wire control for minimizing tower
% top displacement. The agent observes tower displacement and wave elevation.

% Environment Setup
% Define observation and action specifications
obsInfo = rlNumericSpec([3 1], 'LowerLimit', [-inf; -inf; -inf], 'UpperLimit', [inf; inf; inf]);
obsInfo.Name = 'observations';
obsInfo.Description = 'Tower displacement tower velocity and wave elevation';

% Action specification (guy wire control displacements)
actInfo = rlNumericSpec([2 1], 'LowerLimit', [-0.02; -0.02], 'UpperLimit', [0.02; 0.02]);
actInfo.Name = 'guy_wire_controls';
actInfo.Description = 'Control displacements for guy wires 1 and 2 [m]';

% Create custom environment with function handles
env = rlFunctionEnv(obsInfo, actInfo, @MonopileStepFcn, @MonopileResetFcn);

%% DDPG Agent Configuration
% Network parameters
sequenceLength = 100;
hiddenSize = 64;
learningRate = 1e-4;

% % Create critic network with LSTM
% criticNetwork = createCriticLSTMNetwork(obsInfo, actInfo, hiddenSize, sequenceLength);
% critic = rlQValueFunction(criticNetwork, obsInfo, actInfo, ...
%     'ObservationInputNames', 'obs', 'ActionInputNames', 'act');
% 
% % Create actor network with LSTM
% actorNetwork = createActorLSTMNetwork(obsInfo, actInfo, hiddenSize, sequenceLength);
% actor = rlContinuousDeterministicActor(actorNetwork, obsInfo, actInfo, ...
%     'ObservationInputNames', 'obs');

% DDPG agent options
% agentOptions = rlDDPGAgentOptions(...
%     'SampleTime', 0.01, ...
%     'CriticOptimizerOptions', rlOptimizerOptions('LearnRate', 1e-4, 'GradientThreshold', 1), ...
%     'ActorOptimizerOptions', rlOptimizerOptions('LearnRate', 1e-5, 'GradientThreshold', 1), ...
%     'ExperienceBufferLength', 1e6, ...
%     'MiniBatchSize', 64, ...
%     'NumStepsToLookAhead', 1, ...
%     'TargetSmoothFactor', 1e-3, ...
%     'SequenceLength', sequenceLength ...  % Add this line for LSTM networks
% );

% Create DDPG agent
%agent = rlDDPGAgent(actor, critic, agentOptions);

%% Training Options
trainOpts = rlTrainingOptions(...
    'MaxEpisodes', 1000, ...
    'MaxStepsPerEpisode', 200, ...
    'ScoreAveragingWindowLength', 50, ...
    'Verbose', true, ...
    'Plots', 'training-progress', ...
    'StopTrainingCriteria', 'AverageReward', ...
    'StopTrainingValue', -0.01, ...
    'SaveAgentCriteria', 'EpisodeReward', ...
    'SaveAgentValue', -100);

%% Train the Agent
fprintf('Starting DDPG training for monopile control...\n');
trainingStats = train(agent, env, trainOpts);

%% Save Trained Agent
save('trainedDDPGAgent.mat', 'agent', 'trainingStats');
fprintf('Training completed. Agent saved to trainedDDPGAgent.mat\n');

end

function criticNetwork = createCriticLSTMNetwork(obsInfo, actInfo, hiddenSize, sequenceLength)
% Create critic network with LSTM for temporal understanding
% Input dimension is observation + action concatenated
inputDim = obsInfo.Dimension(1) + actInfo.Dimension(1);

criticNetwork = [
    sequenceInputLayer(inputDim, 'Name', 'obs_act')
    lstmLayer(hiddenSize, 'OutputMode', 'sequence', 'Name', 'lstm1')
    dropoutLayer(0.2, 'Name', 'dropout1')
    lstmLayer(32, 'OutputMode', 'sequence', 'Name', 'lstm2')  % Output only last time step
    fullyConnectedLayer(16, 'Name', 'fc1')
    reluLayer('Name', 'relu1')
    fullyConnectedLayer(1, 'Name', 'output')
];
end

function actorNetwork = createActorLSTMNetwork(obsInfo, actInfo, hiddenSize, sequenceLength)
% Create actor network with LSTM for temporal understanding

actorNetwork = [
    sequenceInputLayer(obsInfo.Dimension(1), 'Name', 'obs')
    lstmLayer(hiddenSize, 'OutputMode', 'sequence', 'Name', 'lstm1')
    dropoutLayer(0.2, 'Name', 'dropout1')
    fullyConnectedLayer(32, 'Name', 'fc1')
    reluLayer('Name', 'relu1')
    fullyConnectedLayer(16, 'Name', 'fc2')
    reluLayer('Name', 'relu2')
    fullyConnectedLayer(actInfo.Dimension(1), 'Name', 'output')
    tanhLayer('Name', 'tanh')
    scalingLayer('Name', 'scaling', 'Scale', 0.5)
];
end

    % Nested step function for monopile environment
    function [nextObs, reward, isDone, loggedSignals] = MonopileStepFcn(action, loggedSignals)
        persistent prevAction prevDisplacement actionHistory stepCount
        
        if isempty(prevAction)
            prevAction = [0; 0];
            prevDisplacement = 0;
            actionHistory = zeros(10, 2); % Store last 10 actions
            stepCount = 0;
        end
        
        % Increment step counter
        stepCount = stepCount + 1;
        
        % Extract guy wire controls from action
        deltaL_1 = action(1);
        deltaL_2 = action(2);
        
        % Call monopile simulation function
        [tip_displacement, tip_velocity, tip_acceleration, base_moment, guy_wire_forces] = ...
            monopileSimulinkFunction(deltaL_1, deltaL_2);
        
        % Get current wave elevation (simplified - you may want to implement proper wave tracking)
        wave_elevation = sin(2*pi*stepCount*0.01/10) * 2.5; % Simplified sinusoidal wave
        
        % Prepare observations
        nextObs = [tip_displacement; tip_velocity; wave_elevation];
        
        % Enhanced reward function
        displacement_penalty = -abs(tip_displacement);
        control_effort_penalty = -0.1 * (abs(action(1)) + abs(action(2)));
        
        % Penalize constant actions (encourage dynamic response)
        action_change = abs(action - prevAction);
        stagnation_penalty = -0.05 * exp(-10 * sum(action_change)); % Exponential penalty for no change
        
        % Reward for reducing displacement rate (derivative control)
        displacement_rate = abs(tip_displacement - prevDisplacement);
        rate_reward = 0.1 * max(0, 0.1 - displacement_rate); % Reward for reducing rate
        
        % Penalize oscillatory control (too much action variation)
        actionHistory = [actionHistory(2:end, :); action'];
        action_variance = var(actionHistory);
        oscillation_penalty = -0.02 * sum(action_variance);
        
        reward = displacement_penalty + control_effort_penalty + stagnation_penalty + rate_reward + oscillation_penalty;
        
        % Update history
        prevAction = action;
        prevDisplacement = tip_displacement;

        % Episode termination conditions
        isDone = stepCount >= 2000; % End if too many steps
        
        % Log signals for analysis
        loggedSignals.tip_displacement = tip_displacement;
        loggedSignals.guy_wire_forces = guy_wire_forces;
        loggedSignals.control_action = action;
        loggedSignals.reward = reward;
    end

    % Nested reset function for monopile environment
    function [initialObs, loggedSignals] = MonopileResetFcn()
        % Reset step counter
        stepCount = 0;
        
        % Reset the monopile simulation by calling it with zero inputs
        [tip_displacement, tip_velocity, ~, ~, guy_wire_forces] = monopileSimulinkFunction(0, 0);
        
        % Initial wave elevation
        wave_elevation = 0;
        
        % Initial observations
        initialObs = [tip_displacement; tip_velocity; wave_elevation];
        
        % Initialize logged signals
        loggedSignals = struct();
        loggedSignals.tip_displacement = tip_displacement;
        loggedSignals.guy_wire_forces = guy_wire_forces;
        loggedSignals.control_action = [0; 0];
        loggedSignals.reward = 0;
    end


