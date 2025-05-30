% Create a Simulink model for the monopile simulator with guy wire control
% This version uses an S-Function for fully automated setup

% Model name
modelName = 'MonopileSimulinkModelSFunc';

% Close any existing model with the same name
if bdIsLoaded(modelName)
    close_system(modelName, 0);
end

% Create a new model
sys = new_system(modelName);
open_system(sys);

% Set solver parameters
set_param(modelName, 'Solver', 'ode45', 'StopTime', '60', 'MaxStep', '0.1');

%% Add blocks

% Add S-Function block
sfunctionBlock = add_block('simulink/User-Defined Functions/S-Function', ...
    [modelName, '/Monopile S-Function']);
set_param(sfunctionBlock, 'Position', [300, 200, 400, 300]);
set_param(sfunctionBlock, 'FunctionName', 'monopileSFunction');

% Force Simulink to update the block to recognize ports
set_param(modelName, 'SimulationCommand', 'update');

% Add input blocks for deltaL_1 and deltaL_2
input1 = add_block('simulink/Sources/Constant', [modelName, '/deltaL_1']);
set_param(input1, 'Value', '0');
set_param(input1, 'Position', [50, 180, 100, 210]);

input2 = add_block('simulink/Sources/Constant', [modelName, '/deltaL_2']);
set_param(input2, 'Value', '0');
set_param(input2, 'Position', [50, 280, 100, 310]);

% Add slider gain blocks for interactive control
slider1 = add_block('simulink/Math Operations/Slider Gain', [modelName, '/Slider1']);
set_param(slider1, 'Position', [150, 180, 200, 210]);
set_param(slider1, 'Gain', '0');

slider2 = add_block('simulink/Math Operations/Slider Gain', [modelName, '/Slider2']);
set_param(slider2, 'Position', [150, 280, 200, 310]);
set_param(slider2, 'Gain', '0');

% Add scopes for outputs
scope1 = add_block('simulink/Sinks/Scope', [modelName, '/Tip Displacement']);
set_param(scope1, 'Position', [550, 120, 600, 150]);

scope2 = add_block('simulink/Sinks/Scope', [modelName, '/Tip Velocity']);
set_param(scope2, 'Position', [550, 170, 600, 200]);

scope3 = add_block('simulink/Sinks/Scope', [modelName, '/Tip Acceleration']);
set_param(scope3, 'Position', [550, 220, 600, 250]);

scope4 = add_block('simulink/Sinks/Scope', [modelName, '/Base Moment']);
set_param(scope4, 'Position', [550, 270, 600, 300]);

scope5 = add_block('simulink/Sinks/Scope', [modelName, '/Guy Wire Forces']);
set_param(scope5, 'Position', [550, 320, 600, 350]);

% Add a mux to combine control inputs
muxCtr = add_block('simulink/Signal Routing/Mux', [modelName, '/MuxCtr']);
set_param(muxCtr, 'Inputs', '2');
set_param(muxCtr, 'Position', [300, 200, 320, 280]);

% Add a demux to separate outputs for individual scopes
Demux = add_block('simulink/Signal Routing/Demux', [modelName, '/Demux']);
set_param(Demux, 'Outputs', '5');
set_param(Demux, 'Position', [450, 200, 470, 300]);

% Add a display block
display = add_block('simulink/Sinks/Display', [modelName, '/Display']);
set_param(display, 'Position', [500, 380, 550, 410]);

% Add signal labels
label1 = add_block('simulink/Signal Attributes/Signal Specification', [modelName, '/deltaL1_label']);
set_param(label1, 'Position', [250, 175, 280, 195]);

label2 = add_block('simulink/Signal Attributes/Signal Specification', [modelName, '/deltaL2_label']);
set_param(label2, 'Position', [250, 275, 280, 295]);

%% Connect blocks

% Connect inputs through sliders to mux
add_line(modelName, 'deltaL_1/1', 'Slider1/1');
add_line(modelName, 'Slider1/1', 'deltaL1_label/1');
add_line(modelName, 'deltaL1_label/1', 'MuxCtr/1');

add_line(modelName, 'deltaL_2/1', 'Slider2/1');
add_line(modelName, 'Slider2/1', 'deltaL2_label/1');
add_line(modelName, 'deltaL2_label/1', 'MuxCtr/2');

% Connect mux output to S-Function input
add_line(modelName, 'MuxCtr/1', 'Monopile S-Function/1');

% Connect S-Function output to demux
add_line(modelName, 'Monopile S-Function/1', 'Demux/1');

% Connect demux outputs to scopes
add_line(modelName, 'Demux/1', 'Tip Displacement/1');
add_line(modelName, 'Demux/2', 'Tip Velocity/1');
add_line(modelName, 'Demux/3', 'Tip Acceleration/1');
add_line(modelName, 'Demux/4', 'Base Moment/1');
add_line(modelName, 'Demux/5', 'Guy Wire Forces/1');



% Connect demux to display (showing all outputs)
add_line(modelName, 'Demux/1', 'Display/1');



% Save the model
save_system(modelName);

% Display success message
fprintf('\n=== Simulink Model with S-Function Created Successfully ===\n');
fprintf('Model Name: %s\n\n', modelName);
fprintf('This model is ready to use immediately!\n');
fprintf('\nHow to use:\n');
fprintf('1. Use the sliders to control deltaL_1 and deltaL_2\n');
fprintf('2. Click "Run" to start the simulation\n');
fprintf('3. Monitor outputs in the scope blocks\n\n');
fprintf('Control Ranges:\n');
fprintf('- deltaL_1 (Right Guy Wire): -0.1 to +0.1 meters\n');
fprintf('- deltaL_2 (Left Guy Wire):  -0.1 to +0.1 meters\n');
fprintf('\nPositive values increase tension, negative values decrease tension.\n');
fprintf('\nNote: The first simulation run may take longer as MATLAB compiles the functions.\n');