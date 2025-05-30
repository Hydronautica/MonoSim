% Test script for the monopile Simulink function
% This script tests the function with different control inputs

clc; clear; close all;

% Add path to parent directory to access all functions
addpath('..');

fprintf('Testing Monopile Simulink Function...\n');
fprintf('=====================================\n\n');

%% Test 1: No control input (baseline)
fprintf('Test 1: Baseline (no control)\n');
deltaL_1 = 0;
deltaL_2 = 0;

tic;
[tip_disp_1, tip_vel_1, tip_acc_1, base_moment_1, guy_forces_1] = ...
    monopileSimulinkFunction(deltaL_1, deltaL_2);
time_1 = toc;

fprintf('  deltaL_1 = %.3f m, deltaL_2 = %.3f m\n', deltaL_1, deltaL_2);
fprintf('  Tip displacement: %.4f m\n', tip_disp_1);
fprintf('  Tip velocity: %.4f m/s\n', tip_vel_1);
fprintf('  Tip acceleration: %.4f m/s²\n', tip_acc_1);
fprintf('  Base moment: %.2e N⋅m\n', base_moment_1);
fprintf('  Guy wire forces: %.2e N\n', guy_forces_1);
fprintf('  Computation time: %.2f seconds\n\n', time_1);

%% Test 2: Positive control input (increase tension)
fprintf('Test 2: Increased tension\n');
deltaL_1 = 0.05;  % 5 cm reduction in rest length
deltaL_2 = 0.05;

tic;
[tip_disp_2, tip_vel_2, tip_acc_2, base_moment_2, guy_forces_2] = ...
    monopileSimulinkFunction(deltaL_1, deltaL_2);
time_2 = toc;

fprintf('  deltaL_1 = %.3f m, deltaL_2 = %.3f m\n', deltaL_1, deltaL_2);
fprintf('  Tip displacement: %.4f m\n', tip_disp_2);
fprintf('  Tip velocity: %.4f m/s\n', tip_vel_2);
fprintf('  Tip acceleration: %.4f m/s²\n', tip_acc_2);
fprintf('  Base moment: %.2e N⋅m\n', base_moment_2);
fprintf('  Guy wire forces: %.2e N\n', guy_forces_2);
fprintf('  Computation time: %.2f seconds\n\n', time_2);

%% Test 3: Negative control input (decrease tension)
fprintf('Test 3: Decreased tension\n');
deltaL_1 = -0.02;  % 2 cm increase in rest length
deltaL_2 = -0.02;

tic;
[tip_disp_3, tip_vel_3, tip_acc_3, base_moment_3, guy_forces_3] = ...
    monopileSimulinkFunction(deltaL_1, deltaL_2);
time_3 = toc;

fprintf('  deltaL_1 = %.3f m, deltaL_2 = %.3f m\n', deltaL_1, deltaL_2);
fprintf('  Tip displacement: %.4f m\n', tip_disp_3);
fprintf('  Tip velocity: %.4f m/s\n', tip_vel_3);
fprintf('  Tip acceleration: %.4f m/s²\n', tip_acc_3);
fprintf('  Base moment: %.2e N⋅m\n', base_moment_3);
fprintf('  Guy wire forces: %.2e N\n', guy_forces_3);
fprintf('  Computation time: %.2f seconds\n\n', time_3);

%% Test 4: Asymmetric control (different inputs for each guy wire)
fprintf('Test 4: Asymmetric control\n');
deltaL_1 = 0.03;   % Right guy wire: increase tension
deltaL_2 = -0.01;  % Left guy wire: decrease tension

tic;
[tip_disp_4, tip_vel_4, tip_acc_4, base_moment_4, guy_forces_4] = ...
    monopileSimulinkFunction(deltaL_1, deltaL_2);
time_4 = toc;

fprintf('  deltaL_1 = %.3f m, deltaL_2 = %.3f m\n', deltaL_1, deltaL_2);
fprintf('  Tip displacement: %.4f m\n', tip_disp_4);
fprintf('  Tip velocity: %.4f m/s\n', tip_vel_4);
fprintf('  Tip acceleration: %.4f m/s²\n', tip_acc_4);
fprintf('  Base moment: %.2e N⋅m\n', base_moment_4);
fprintf('  Guy wire forces: %.2e N\n', guy_forces_4);
fprintf('  Computation time: %.2f seconds\n\n', time_4);

%% Summary and Analysis
fprintf('Summary and Analysis\n');
fprintf('===================\n');
fprintf('Control Effect on Tip Displacement:\n');
fprintf('  Baseline:           %.4f m\n', tip_disp_1);
fprintf('  Increased tension:  %.4f m (%.1f%% change)\n', tip_disp_2, 100*(tip_disp_2-tip_disp_1)/abs(tip_disp_1));
fprintf('  Decreased tension:  %.4f m (%.1f%% change)\n', tip_disp_3, 100*(tip_disp_3-tip_disp_1)/abs(tip_disp_1));
fprintf('  Asymmetric control: %.4f m (%.1f%% change)\n', tip_disp_4, 100*(tip_disp_4-tip_disp_1)/abs(tip_disp_1));

fprintf('\nControl Effect on Guy Wire Forces:\n');
fprintf('  Baseline:           %.2e N\n', guy_forces_1);
fprintf('  Increased tension:  %.2e N (%.1f%% change)\n', guy_forces_2, 100*(guy_forces_2-guy_forces_1)/guy_forces_1);
fprintf('  Decreased tension:  %.2e N (%.1f%% change)\n', guy_forces_3, 100*(guy_forces_3-guy_forces_1)/guy_forces_1);
fprintf('  Asymmetric control: %.2e N (%.1f%% change)\n', guy_forces_4, 100*(guy_forces_4-guy_forces_1)/guy_forces_1);

fprintf('\nAverage computation time: %.2f seconds\n', mean([time_1, time_2, time_3, time_4]));

%% Create comparison plot
figure('Position', [100, 100, 1200, 800]);

% Tip displacement comparison
subplot(2, 3, 1);
bar([tip_disp_1, tip_disp_2, tip_disp_3, tip_disp_4]);
set(gca, 'XTickLabel', {'Baseline', 'Inc. Tension', 'Dec. Tension', 'Asymmetric'});
ylabel('Tip Displacement [m]');
title('Control Effect on Tip Displacement');
grid on;

% Tip velocity comparison
subplot(2, 3, 2);
bar([tip_vel_1, tip_vel_2, tip_vel_3, tip_vel_4]);
set(gca, 'XTickLabel', {'Baseline', 'Inc. Tension', 'Dec. Tension', 'Asymmetric'});
ylabel('Tip Velocity [m/s]');
title('Control Effect on Tip Velocity');
grid on;

% Tip acceleration comparison
subplot(2, 3, 3);
bar([tip_acc_1, tip_acc_2, tip_acc_3, tip_acc_4]);
set(gca, 'XTickLabel', {'Baseline', 'Inc. Tension', 'Dec. Tension', 'Asymmetric'});
ylabel('Tip Acceleration [m/s²]');
title('Control Effect on Tip Acceleration');
grid on;

% Base moment comparison
subplot(2, 3, 4);
bar([base_moment_1, base_moment_2, base_moment_3, base_moment_4]);
set(gca, 'XTickLabel', {'Baseline', 'Inc. Tension', 'Dec. Tension', 'Asymmetric'});
ylabel('Base Moment [N⋅m]');
title('Control Effect on Base Moment');
grid on;

% Guy wire forces comparison
subplot(2, 3, 5);
bar([guy_forces_1, guy_forces_2, guy_forces_3, guy_forces_4]);
set(gca, 'XTickLabel', {'Baseline', 'Inc. Tension', 'Dec. Tension', 'Asymmetric'});
ylabel('Guy Wire Forces [N]');
title('Control Effect on Guy Wire Forces');
grid on;

% Control inputs visualization
subplot(2, 3, 6);
control_inputs = [0, 0; 0.05, 0.05; -0.02, -0.02; 0.03, -0.01];
bar(control_inputs);
set(gca, 'XTickLabel', {'Baseline', 'Inc. Tension', 'Dec. Tension', 'Asymmetric'});
ylabel('Control Input [m]');
title('Control Inputs (deltaL)');
legend('deltaL_1 (Right)', 'deltaL_2 (Left)', 'Location', 'best');
grid on;

sgtitle('Monopile Simulink Function Test Results', 'FontSize', 16, 'FontWeight', 'bold');

fprintf('\nTest completed successfully!\n');
fprintf('Plots generated showing control effectiveness.\n');