function plotSurfaceElevation(time, eta2, params)
% PLOTSURFACEELEVATION Plots wave surface elevation
%
% Inputs:
%   time   - Time vector
%   eta2   - Second-order surface elevation component
%   params - Structure containing simulation parameters

% Extract parameters
Nfreq = params.Nfreq;
omegaVec = params.omegaVec;
A_m = params.A_m;
phiRand = params.phiRand;
% Ensure vectors have compatible dimensions for element-wise operations
if size(omegaVec, 1) > 1 && size(omegaVec, 2) == 1
    % Convert column vector to row vector if needed
    omegaVec = omegaVec';
end

if size(A_m, 1) > 1 && size(A_m, 2) == 1
    % Convert column vector to row vector if needed
    A_m = A_m';
end

if size(phiRand, 1) > 1 && size(phiRand, 2) == 1
    % Convert column vector to row vector if needed
    phiRand = phiRand';
end
% Calculate first-order irregular wave elevation
eta_irreg = zeros(1, length(time));
for i = 1:length(time)
    t = time(i);
    eta_irreg(i) = sum(A_m .* cos(omegaVec*t + phiRand));
end

% Plot both first-order and total (first + second order) surface elevation
figure; hold on; grid on;
plot(time, eta_irreg, 'LineWidth', 2, 'DisplayName', 'First-order');

% Ensure eta2 is the same size as eta_irreg for plotting
if length(eta2) == 1 && length(eta_irreg) > 1
    % If eta2 is a scalar, expand it to match eta_irreg
    eta2_plot = ones(size(eta_irreg)) * eta2;
else
    % Otherwise use as is, assuming it's already the right size
    eta2_plot = eta2;
end

plot(time, eta_irreg + eta2_plot, 'LineWidth', 2, 'DisplayName', 'Total (1st + 2nd order)');
xlabel('Time [s]');
ylabel('Surface elevation \eta [m]');
title('Wave surface elevation time-series');
legend('Location', 'Best');

end