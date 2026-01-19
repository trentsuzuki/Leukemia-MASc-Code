%% Generate plots for AMPC + UKF for thesis

clear; clc; close all;

set(0,'defaultTextInterpreter','latex');

% Load results from a previous script run
% filename = "ga_ukf_fullJostModel_20250811_0148";
% filename = "ga_ukf_3param_est_20250818_0410";
% filename = "ga_ukf_2param_est_20250819_0943";
filename = "ga_ukf_fullJostModel_20250827_1021";
load("results/" + filename + ".mat");

%% Parameters

y_r = ancTarget; % Desired steady-state neutrophil concentration (G/L)
neutrophil_color = [0.5 0.5 0.9 0.2];
drug_effect_color = [0.5 0.5 0.9 0.2];
y_r_color = [0.9 0.14 0.14];

%% Generate plots
fig_traj = figure('Visible', 'on');
set(fig_traj, 'Color', 'White', 'Theme', 'Light');

subplot(2,1,1);
hold on;
title('AMPC: Simulated Treatment Outcomes', ...
    'FontSize', 16);
xlabel('Time (days)');
ylabel('ANC (G/L)');
grid on;
xlim([0 max(allResults{1}.t)])
tmp = cellfun(@(p) p.x(:,nx), allResults, 'UniformOutput', false);
allAncMeas = vertcat(tmp{:});
ylim([0 1.3*max(allAncMeas)])

subplot(2,1,2);
hold on;
xlabel('Time (days)');
ylabel('Daily 6-MP Dose (mg)');
grid on;
xlim([0 max(allResults{1}.t)])
ylim([0 200])

for i = 1:length(allResults)
    id = allResults{i}.patientID; % Extract patient ID
    t = allResults{i}.t; % Extract time vector for the current patient
    y = allResults{i}.x(:,nx); % Extract neutrophil measurements
    u = allResults{i}.u; % Extract drug doses

    % Subplot for Neutrophil Concentration
    subplot(2,1,1);
    % title(['Patient ', num2str(id), ' Treatment Outcome'], 'FontSize', 16);
    plot(t, y, ...
        '-', ...
        'Color', neutrophil_color, ...
        'LineWidth', 1, ...
        'HandleVisibility', 'off');

    % Subplot for Drug Effect
    subplot(2,1,2);
    t_zoh = [t(1); kron(t(2:end)', [1; 1])];
    u_zoh = [kron(u(1:end-1), [1; 1]); u(end-1);];
    d = plot(t_zoh, u_zoh, '-', ...
        'Color', drug_effect_color, ...
        'Linewidth', 2, ...
        'HandleVisibility', 'off');
end

subplot(2,1,1)
ref = yline(y_r, '--', 'Color', y_r_color, 'LineWidth', 2, 'DisplayName', 'Reference (y_r)');
legend('Interpreter', 'latex', 'TextColor', y_r_color)
trajectory_filename = sprintf("figs/" + filename + "_allPatients_Outcome.png");
exportgraphics(fig_traj, trajectory_filename, 'Resolution', 300);
close(fig_traj)

