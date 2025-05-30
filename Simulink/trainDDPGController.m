function trainDDPGController()
% TRAINDDPGCONTROLLER Train DDPG agent with LSTM networks for monopile control
%
% This function implements a Deep Deterministic Policy Gradient (DDPG) agent
% with LSTM networks to learn optimal guy wire control for minimizing tower
% top displacement. The agent observes tower displacement and wave elevation.

%% Environment Setup
% Define observation and action specifications
obsInfo = rlNumericSpec([2 1], 'LowerLimit', [-inf; -inf], 'UpperLimit', [inf; inf]);
obsInfo.Name = 'observations';
obsInfo.Description = 'Tower displacement and wave elevation';

% Action specification (guy wire control displacements)
actInfo = rlNumericSpec([2 1], 'LowerLimit', [-0.002; -0.002], 'UpperLimit', [0.002; 0.002]);
actInfo.Name = 'guy_wire_controls';
actInfo.Description = 'Control displacements for guy wires 1 and 2 [m]';

% Create custom environment with function handles
env = rlFunctionEnv(obsInfo, actInfo, @MonopileStepFcn, @MonopileResetFcn);

%% DDPG Agent Configuration
% Network parameters
sequenceLength = 100;
hiddenSize = 64;
learningRate = 1e-4;

% Create critic network with LSTM
criticNetwork = createCriticLSTMNetwork(obsInfo, actInfo, hiddenSize, sequenceLength);
critic = rlQValueFunction(criticNetwork, obsInfo, actInfo, ...
    'ObservationInputNames', 'obs', 'ActionInputNames', 'act');

% Create actor network with LSTM
actorNetwork = createActorLSTMNetwork(obsInfo, actInfo, hiddenSize, sequenceLength);
actor = rlContinuousDeterministicActor(actorNetwork, obsInfo, actInfo, ...
    'ObservationInputNames', 'obs');

% DDPG agent options
agentOptions = rlDDPGAgentOptions(...
    'SampleTime', 0.01, ...
    'CriticOptimizerOptions', rlOptimizerOptions('LearnRate', learningRate), ...
    'ActorOptimizerOptions', rlOptimizerOptions('LearnRate', learningRate/10), ...
    'ExperienceBufferLength', 1e6, ...
    'MiniBatchSize', 128, ...
    'NumStepsToLookAhead', 10, ...
    'TargetSmoothFactor', 1e-3 ...
    );

% Create DDPG agent
agent = rlDDPGAgent(actor, critic, agentOptions);

%% Training Options
trainOpts = rlTrainingOptions(...
    'MaxEpisodes', 100, ...
    'MaxStepsPerEpisode', 2000, ...
    'ScoreAveragingWindowLength', 50, ...
    'Verbose', true, ...
    'Plots', 'training-progress', ...
    'StopTrainingCriteria', 'AverageReward', ...
    'StopTrainingValue', -0.1, ...
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
% Create critic network with MLP layers (no LSTM)

% Observation path
obsPath = [
    featureInputLayer(obsInfo.Dimension(1), 'Name', 'obs')
    fullyConnectedLayer(hiddenSize, 'Name', 'fc_obs1')
    reluLayer('Name', 'relu_obs1')
    fullyConnectedLayer(32, 'Name', 'fc_obs2')
    reluLayer('Name', 'relu_obs2')
];

% Action path
actionPath = [
    featureInputLayer(actInfo.Dimension(1), 'Name', 'act')
    fullyConnectedLayer(32, 'Name', 'fc_act')
    reluLayer('Name', 'relu_act')
];

% Common path (combine observation and action)
commonPath = [
    additionLayer(2, 'Name', 'add')
    fullyConnectedLayer(64, 'Name', 'fc_common1')
    reluLayer('Name', 'relu_common1')
    fullyConnectedLayer(32, 'Name', 'fc_common2')
    reluLayer('Name', 'relu_common2')
    fullyConnectedLayer(1, 'Name', 'output')
];

% Create layer graph
criticNetwork = layerGraph(obsPath);
criticNetwork = addLayers(criticNetwork, actionPath);
criticNetwork = addLayers(criticNetwork, commonPath);

% Connect layers
criticNetwork = connectLayers(criticNetwork, 'relu_obs2', 'add/in1');
criticNetwork = connectLayers(criticNetwork, 'relu_act', 'add/in2');

end

function actorNetwork = createActorLSTMNetwork(obsInfo, actInfo, hiddenSize, sequenceLength)
% Create actor network with MLP layers (no LSTM)

actorNetwork = [
    featureInputLayer(obsInfo.Dimension(1), 'Name', 'obs')
    fullyConnectedLayer(hiddenSize, 'Name', 'fc1')
    reluLayer('Name', 'relu1')
    fullyConnectedLayer(32, 'Name', 'fc2')
    reluLayer('Name', 'relu2')
    fullyConnectedLayer(16, 'Name', 'fc3')
    reluLayer('Name', 'relu3')
    fullyConnectedLayer(actInfo.Dimension(1), 'Name', 'output')
    tanhLayer('Name', 'tanh')
    scalingLayer('Name', 'scaling', 'Scale', 0.5) % Scale to [-0.5, 0.5] range
];

end

    % Nested step function for monopile environment
    function [nextObs, reward, isDone, loggedSignals] = MonopileStepFcn(action, loggedSignals)
        persistent simData stepCount
        
        if isempty(stepCount)
            stepCount = 0;
        end
        
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
        nextObs = [tip_displacement; wave_elevation];
        
        % Calculate reward (negative displacement to minimize it)
        %reward = -abs(tip_displacement) - 0.1*(abs(deltaL_1) + abs(deltaL_2)); % Penalize large control actions
        reward = -abs(tip_displacement)^2; % Penalize large control actions

        % Episode termination conditions
        isDone = stepCount >= 2000 ; % End if too many steps or excessive displacement
        
        % Log signals for analysis
        loggedSignals.tip_displacement = tip_displacement;
        loggedSignals.guy_wire_forces = guy_wire_forces;
        loggedSignals.control_action = action;
        loggedSignals.reward = reward;
    end

    % Nested reset function for monopile environment
    function [initialObs, loggedSignals] = MonopileResetFcn()
        persistent stepCount
        stepCount = 0;
        
        % Reset the monopile simulation by calling it with zero inputs
        [tip_displacement, ~, ~, ~, guy_wire_forces] = monopileSimulinkFunction(0, 0);
        
        % Initial wave elevation
        wave_elevation = 0;
        
        % Initial observations
        initialObs = [tip_displacement; wave_elevation];
        
        % Initialize logged signals
        loggedSignals = struct();
        loggedSignals.tip_displacement = tip_displacement;
        loggedSignals.guy_wire_forces = guy_wire_forces;
        loggedSignals.control_action = [0; 0];
        loggedSignals.reward = 0;
    end

