% This script runs Monte Carlo simulations for the Friberg model for multiple
% patients. It reads patient-specific parameters (Base, ktr, gamma) from a
% CSV file and runs 1000 simulations for each patient.

clear; clc; close all;
set(0,'defaultTextInterpreter','latex');

% --- Setup Paths ---
userdir = string(java.lang.System.getProperty('user.home')) + "/";
% NOTE: The following paths should be updated to match your system.
codeFolder = userdir + "UofT/Research/Leukemia/Code/";
% Path to CSV file with patient params (e.g., Data Sheet 1 from Jost 2020 supplementary materials)
paramFilePath = codeFolder + "Jost 2020 Data Analysis/Supplementary Materials/Data Sheet 1.CSV";

%% Read Patient Data from CSV
% The script expects a CSV file with columns: ID, Base, ktr, gamma
try
    patientData = readtable(paramFilePath);
    fprintf('Successfully loaded %d patients from CSV\n\n', height(patientData));
catch ME
    if strcmp(ME.identifier, 'MATLAB:readtable:FileNotFound')
        error('Error: Patient data CSV file not found at the specified path. Please check paramFilePath.');
    else
        rethrow(ME);
    end
end

%% Parameters
% --- Fixed Model & Control Parameters ---
k_ma = 2.375;      % Maturation rate constant (1/day)
y_r = 1.5;         % Desired steady-state neutrophil concentration (G/L)
y_min_safe = 0.5;  % Safe lower limit in clinical practice
eps = 2.1 * y_r;   % Tuning parameter for control law

% --- Simulation Parameters ---
num_simulations = 100; % Number of simulations to run per patient
t_span = [0 300];       % Time span for each simulation (days)
sd = 0.75;              % Standard deviation for randomizing initial conditions

% --- Plot options ---
neutrophil_color = [0.5 0.5 0.9 0.01]; % Blue with transparency
drug_effect_color = [0.5 0.5 0.9 0.01]; % Blue with transparency
y_r_color = [0.9 0.14 0.14];         % Red for reference line

%% IMPORTANT: Options to save plots as files
save_trajectory_figures = true;

%% Prepare Figure for Combined Trajectory Plot
% Create the figure and subplots *before* looping through patients.
fig_traj = figure('Visible', 'off'); % Set to 'on' to see it build, 'off' for background processing
set(fig_traj, 'Color', 'White', 'Theme', 'Light', 'Position', [100 100 900 700]);

% Subplot for Neutrophil Concentration
subplot(2,1,1);
hold on;
title(['Output Feedback Control: All Patients  (', ...
    num2str(num_simulations), ' simulations per patient)'], 'FontSize', 16);
xlabel('Time (days)');
ylabel('Neutrophil Conc. (y)');
grid on;

% Subplot for Drug Effect
subplot(2,1,2);
hold on;
xlabel('Time (days)');
ylabel('Drug Effect (u)');
grid on;

% To track the overall maximum y-value for final plot scaling
overall_max_y = 0;

%% Main Loop for Each Patient
for p_idx = 1:height(patientData) 
    % --- Get current patient's data ---
    current_patient = patientData(p_idx, :);
    patient_id = current_patient.ID;
    
    fprintf('--- Starting Simulations for Patient %d ---\n', patient_id);
    fprintf('Parameters: Base=%.3f, ktr=%.3f, gamma=%.3f\n', ...
        current_patient.Base, current_patient.ktr, current_patient.gamma);
        
    %% Model & Simulation Parameters
    % --- Model Parameters (loaded from file) ---
    k_tr = current_patient.ktr;   % Transit rate constant (1/day)
    Base = current_patient.Base;  % Baseline neutrophil concentration (G/L)
    gamma = current_patient.gamma;% Exponent for homeostatic response
    
    % Store parameters in a struct to pass to the ODE function
    params.k_tr = k_tr;
    params.k_ma = k_ma;
    params.Base = Base;
    params.gamma = gamma;
    
    % Calculate patient-specific equilibrium
    X_r = (k_ma / k_tr) * y_r; % Corresponding equilibrium for precursor cells
    
    %% Monte Carlo Simulation Loop
    % Pre-allocate storage for this patient's results
    final_y_values = zeros(num_simulations, 1);
    max_y_values = zeros(num_simulations, 1);
    min_y_values = zeros(num_simulations, 1);
    
    % Define the base initial state vector for the current patient
    X0 = (k_ma / k_tr) * Base;
    x0_base = [X0; X0; X0; X0; Base];
    
    % --- Inner Simulation Loop ---
    for i = 1:num_simulations
        % Generate random initial conditions (must be positive)
        x0 = -1;
        while any(x0 <= 0)
            x0 = x0_base .* (1 + sd * randn(5, 1));
        end
        
        % Define the (saturated) output feedback control law
        u = @(x) max(0, 1 - (y_r + eps)./(x(5,:) + eps) .* (x(5,:)./Base).^gamma);
        
        % Simulate the system using ode45
        [t, x] = ode45(@(t,x) friberg_model(t, x, u(x), params), t_span, x0);
        
        % Extract the output (neutrophil concentration)
        y = x(:, 5);
        
        % Plot the output trajectory on the shared figure
        subplot(2,1,1);
        plot(t, y, 'Color', neutrophil_color, 'HandleVisibility', 'off');
        
        % Plot the control input on the shared figure
        subplot(2,1,2);
        plot(t, u(x'), 'Color', drug_effect_color, 'HandleVisibility', 'off');
        
        % Store results for this patient's statistics
        final_y_values(i) = y(end);
        max_y_values(i) = max(y);
        min_y_values(i) = min(y);
    end
    
    % Update the overall maximum y-value seen so far
    overall_max_y = max(overall_max_y, max(max_y_values));
end
fprintf('--- All patients processed. ---\n');

%% --- Finalize and Save Combined Trajectory Plot ---
% Add reference line and adjust axis now that all data is plotted
subplot(2,1,1);
yline(y_r, '--', 'Color', y_r_color, 'LineWidth', 2, 'DisplayName', 'Desired Output (y_r)');
text(t_span(2)*0.9, y_r + 0.5, 'Reference', 'Color', y_r_color, 'FontSize', 12);
ylim([0, overall_max_y * 1.1]); % Adjust y-axis based on all simulations

% Save the combined figure
if save_trajectory_figures
    trajectory_filename = 'figs/Output_fb_All_Patients_Combined_Trajectories.png';
    exportgraphics(fig_traj, trajectory_filename, 'Resolution', 300);
    fprintf('Saved combined trajectory plot to: %s\n', trajectory_filename);
end

% close(fig_traj);
