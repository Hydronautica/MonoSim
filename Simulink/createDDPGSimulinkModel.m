function createDDPGSimulinkModel()
% CREATEDDPGSIMULINKMODEL Create a Simulink model with DDPG controller
%
% This function creates a Simulink model that integrates the trained DDPG
% controller with the monopile simulation for real-time control using S-Functions.

%% Check for trained agent
try
    load('trainedDDPGAgent.mat', 'agent');
    fprintf('Loaded trained DDPG agent successfully.\n');
catch
    warning('Could not load trained agent. Creating model with placeholder for agent.');
    agent = [];
end

%% Create new Simulink model
modelName = 'MonopileDDPGControl';

% Close model if it exists
if bdIsLoaded(modelName)
    close_system(modelName, 0);
end

% Create new model
sys = new_system(modelName);
open_system(sys);

% Set solver parameters
set_param(modelName, 'Solver', 'ode45', 'StopTime', '200', 'MaxStep', '0.1');

%% Add blocks

% Add S-Function block for monopile simulation
sfunctionBlock = add_block('simulink/User-Defined Functions/S-Function', ...
    [modelName, '/Monopile_SFunction']);
set_param(sfunctionBlock, 'Position', [400, 200, 500, 300]);
set_param(sfunctionBlock, 'FunctionName', 'monopileSFunction');

% Force Simulink to update the block to recognize ports
set_param(modelName, 'SimulationCommand', 'update');

% Wave generator
waveGen = add_block('simulink/Sources/Sine Wave', [modelName, '/Wave_Generator']);
set_param(waveGen, 'Amplitude', '2.5', 'Frequency', '0.1', 'Phase', '0');
set_param(waveGen, 'Position', [50, 100, 100, 130]);

% DDPG Controller S-Function (we'll create this)
ddpgBlock = add_block('simulink/User-Defined Functions/S-Function', ...
    [modelName, '/DDPG_Controller']);
set_param(ddpgBlock, 'Position', [250, 200, 350, 280]);
set_param(ddpgBlock, 'FunctionName', 'ddpgControllerSFunction');

% Reset signal
resetBlock = add_block('simulink/Sources/Constant', [modelName, '/Reset']);
set_param(resetBlock, 'Value', '0', 'Position', [50, 300, 100, 330]);

% Add Mux for controller inputs (displacement, wave, reset)
muxController = add_block('simulink/Signal Routing/Mux', [modelName, '/Mux_Controller']);
set_param(muxController, 'Inputs', '3', 'Position', [150, 200, 170, 280]);

% Add Mux for monopile inputs (deltaL_1, deltaL_2)
muxMonopile = add_block('simulink/Signal Routing/Mux', [modelName, '/Mux_Monopile']);
set_param(muxMonopile, 'Inputs', '2', 'Position', [350, 220, 370, 260]);

% Add Demux for monopile outputs
demuxMonopile = add_block('simulink/Signal Routing/Demux', [modelName, '/Demux_Monopile']);
set_param(demuxMonopile, 'Outputs', '5', 'Position', [550, 200, 570, 300]);

% Add scopes for visualization
scope1 = add_block('simulink/Sinks/Scope', [modelName, '/Tip_Displacement']);
set_param(scope1, 'Position', [650, 120, 700, 150]);

scope2 = add_block('simulink/Sinks/Scope', [modelName, '/Tip_Velocity']);
set_param(scope2, 'Position', [650, 170, 700, 200]);

scope3 = add_block('simulink/Sinks/Scope', [modelName, '/Tip_Acceleration']);
set_param(scope3, 'Position', [650, 220, 700, 250]);

scope4 = add_block('simulink/Sinks/Scope', [modelName, '/Base_Moment']);
set_param(scope4, 'Position', [650, 270, 700, 300]);

scope5 = add_block('simulink/Sinks/Scope', [modelName, '/Guy_Wire_Forces']);
set_param(scope5, 'Position', [650, 320, 700, 350]);

% Add control scopes
controlScope = add_block('simulink/Sinks/Scope', [modelName, '/Control_Actions']);
set_param(controlScope, 'NumInputPorts', '2', 'Position', [650, 380, 700, 420]);

% Add To Workspace blocks for data logging
logDisp = add_block('simulink/Sinks/To Workspace', [modelName, '/Log_Displacement']);
set_param(logDisp, 'VariableName', 'displacement_data', 'Position', [750, 120, 800, 150]);

logControl = add_block('simulink/Sinks/To Workspace', [modelName, '/Log_Control']);
set_param(logControl, 'VariableName', 'control_data', 'Position', [750, 380, 800, 410]);

logWave = add_block('simulink/Sinks/To Workspace', [modelName, '/Log_Wave']);
set_param(logWave, 'VariableName', 'wave_data', 'Position', [750, 50, 800, 80]);

%% Connect blocks

% Connect inputs to controller mux
add_line(modelName, 'Wave_Generator/1', 'Mux_Controller/2');
add_line(modelName, 'Reset/1', 'Mux_Controller/3');

% Connect controller mux to DDPG controller
add_line(modelName, 'Mux_Controller/1', 'DDPG_Controller/1');

% Connect DDPG outputs to monopile mux
add_line(modelName, 'DDPG_Controller/1', 'Mux_Monopile/1');
add_line(modelName, 'DDPG_Controller/2', 'Mux_Monopile/2');

% Connect monopile mux to S-Function
add_line(modelName, 'Mux_Monopile/1', 'Monopile_SFunction/1');

% Connect S-Function output to demux
add_line(modelName, 'Monopile_SFunction/1', 'Demux_Monopile/1');

% Feedback connection: tip displacement to controller
add_line(modelName, 'Demux_Monopile/1', 'Mux_Controller/1');

% Connect demux outputs to scopes
add_line(modelName, 'Demux_Monopile/1', 'Tip_Displacement/1');
add_line(modelName, 'Demux_Monopile/2', 'Tip_Velocity/1');
add_line(modelName, 'Demux_Monopile/3', 'Tip_Acceleration/1');
add_line(modelName, 'Demux_Monopile/4', 'Base_Moment/1');
add_line(modelName, 'Demux_Monopile/5', 'Guy_Wire_Forces/1');

% Connect control actions to scope
add_line(modelName, 'DDPG_Controller/1', 'Control_Actions/1');
add_line(modelName, 'DDPG_Controller/2', 'Control_Actions/2');

% Connect to workspace blocks
add_line(modelName, 'Demux_Monopile/1', 'Log_Displacement/1');
add_line(modelName, 'Mux_Monopile/1', 'Log_Control/1');
add_line(modelName, 'Wave_Generator/1', 'Log_Wave/1');

%% Force update and arrange
set_param(modelName, 'SimulationCommand', 'update');

%% Save model
save_system(modelName);
fprintf('Created DDPG Simulink model: %s\n', modelName);
fprintf('Note: You need to create ddpgControllerSFunction.m for the DDPG controller\n');

end