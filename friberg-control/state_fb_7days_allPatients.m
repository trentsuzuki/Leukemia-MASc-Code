% This script runs Monte Carlo simulations for the Friberg model for multiple
% patients using a 7-day hold state feedback control method. It reads 
% patient-specific parameters (Base, ktr, gamma) from a CSV file and runs 
% multiple simulations for each patient. All trajectories are overlaid onto 
% a single figure to visualize the overall system performance.

clear; clc; close all;
set(0,'defaultTextInterpreter','latex');

% --- Setup Paths ---
userdir = string(java.lang.System.getProperty('user.home')) + "/";
% NOTE: The following path should be updated to match your system.
codeFolder = userdir + "UofT/Research/Leukemia/Code/";
% Path to CSV file with patient params (e.g., Data Sheet 1 from Jost 2020 supplementary materials)
paramFilePath = codeFolder + "Jost 2020 Data Analysis/Supplementary Materials/Data Sheet 1.CSV";

% NOTE: This script requires the 'friberg_model.m' function to be in the 
% MATLAB path. The function should have the signature:
%       dxdt = friberg_model(t, x, u, params)
% where 'params' is a struct containing k_tr, k_ma, Base, and gamma.

%% Read Patient Data from CSV
% The script expects a CSV file with columns: ID, Base, ktr, gamma
try
    patientData = readtable(paramFilePath);
    fprintf('Successfully loaded %d patients from CSV\n\n', height(patientData));
catch ME
    if strcmp(ME.identifier, 'MATLAB:readtable:FileNotFound')
        error('Error: Patient data CSV file not found. Please check the path in `paramFilePath`.');
    else
        rethrow(ME);
    end
end

%% Parameters
% --- Fixed Model & Control Parameters ---
k_ma = 2.375;      % Maturation rate constant (1/day)
y_r = 1.5;         % Desired steady-state neutrophil concentration (G/L)

% --- Simulation Parameters ---
num_simulations = 50; % Number of simulations to run per patient
t_span = [0 200];       % Time span for each simulation (days)
sd = 0.75;              % Standard deviation for randomizing initial conditions
measurement_interval = 7; % Interval for control update (days)

% --- Plot options ---
neutrophil_color = [0.5 0.5 0.9 0.02];  % Blue with high transparency
drug_effect_color = [0.5 0.5 0.9 0.02]; % Blue with high transparency
y_r_color = [0.9 0.14 0.14];            % Red for reference line
ANCsafeUpperLim = 10;                   % Max safe ANC

%% IMPORTANT: Option to save the final plot
save_trajectory_figure = true;

%% Prepare Figure for Combined Trajectory Plot
% Create the figure and subplots *before* looping through patients.
fig_traj = figure('Visible', 'off'); % Set to 'on' to see it build, 'off' for background processing
set(fig_traj, 'Color', 'White', 'Theme', 'Light', 'Position', [100 100 900 700]);

% Subplot for Neutrophil Concentration
subplot(2,1,1);
hold on;
title(['7-Day Hold State Feedback Control: All Patients (', ...
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
        
    %% Model Parameters for Current Patient
    k_tr = current_patient.ktr;   % Transit rate constant (1/day)
    Base = current_patient.Base;  % Baseline neutrophil concentration (G/L)
    gamma = current_patient.gamma;% Exponent for homeostatic response
    
    % Store parameters in a struct to pass to the ODE function
    params.k_tr = k_tr;
    params.k_ma = k_ma;
    params.Base = Base;
    params.gamma = gamma;
    
    % Define the base initial state vector for the current patient
    X0 = (k_ma / k_tr) * Base;
    x0_base = [X0; X0; X0; X0; Base];
    
    %% Monte Carlo Simulation Loop (Inner Loop)
    for i = 1:num_simulations
        % Generate random initial conditions (must be positive)
        x0 = -1;
        while any(x0 <= 0)
            x0 = x0_base .* (1 + sd * randn(5, 1));
        end
        
        % --- 7-Day Hold Simulation Logic ---
        t_current = t_span(1);
        x_current = x0;
        
        T_full = []; % Stores the complete time vector for this simulation
        X_full = []; % Stores the complete state trajectory
        U_full = []; % Stores the control input for plotting
        
        % Define the control law. It takes a single state vector as input.
        u_law = @(x_vec) max(0, 1 - (x_vec(4) * y_r) / (x_vec(5) * x_vec(1)) * (x_vec(5) / Base)^gamma);
        
        while t_current < t_span(2)
            % Calculate control input based on the state at the start of the interval
            u_val = u_law(x_current);
            
            % Define the time span for the next simulation segment
            t_end = min(t_current + measurement_interval, t_span(2));
            t_interval = [t_current, t_end];
            
            % Simulate the system for this segment using the fixed control input
            [t_segment, x_segment] = ode45(@(t, x) friberg_model(t, x, u_val, params), t_interval, x_current);
            
            % Aggregate results, avoiding duplicate time points at junctions
            if isempty(T_full)
                T_full = t_segment;
                X_full = x_segment;
                U_full = ones(size(t_segment)) * u_val;
            else
                % Append from the second element to avoid duplicating the junction point
                T_full = [T_full; t_segment(2:end)];
                X_full = [X_full; x_segment(2:end, :)];
                U_full = [U_full; ones(length(t_segment) - 1, 1) * u_val];
            end
            
            % Update the state and time for the next iteration
            t_current = t_end;
            x_current = x_segment(end, :)'; % Get the last state as a column vector
        end
        
        y = X_full(:, 5); % Extract neutrophil concentration
        
        % --- Plotting on the Combined Figure ---
        % Plot the neutrophil trajectory
        subplot(2,1,1);
        if max(y) > ANCsafeUpperLim
            plot_color = [neutrophil_color(1:3) 0.25];
        else
            plot_color = neutrophil_color;
        end
        plot(T_full, y, 'Color', plot_color, 'HandleVisibility', 'off');
        
        % Plot the piecewise-constant control input
        subplot(2,1,2);
        % Create points for a stair-step plot
        t_zoh = [T_full(1); kron(T_full(2:end), [1; 1])];
        u_zoh = [kron(U_full(1:end-1), [1; 1]); U_full(end-1);];
        plot(t_zoh, u_zoh, '-', 'Color', drug_effect_color, 'HandleVisibility', 'off');

        % Update the overall maximum y-value seen so far for final scaling
        overall_max_y = max(overall_max_y, max(y));
    end
end
fprintf('--- All patients processed. ---\n');

%% --- Finalize and Save Combined Trajectory Plot ---
% Add reference line and adjust axes now that all data is plotted
subplot(2,1,1);
yline(y_r, '--', 'Color', y_r_color, 'LineWidth', 2, 'DisplayName', 'Desired Output ($y_r$)');
legend('show', 'Location', 'northeast', 'Interpreter', 'latex', 'TextColor', y_r_color);
% text(t_span(2)*0.85, y_r + 0.5, 'Reference', 'Color', y_r_color, 'FontSize', 12);
ylim([0, overall_max_y * 1.1]); % Adjust y-axis based on all simulations

% Set y-limit for the control input plot
subplot(2,1,2);
ylim([0 1.05]);

% Save the combined figure
if save_trajectory_figure
    % Create 'figs' directory if it doesn't exist
    if ~exist('figs', 'dir')
       mkdir('figs')
    end
    trajectory_filename = 'figs/7Day_Hold_State_fb_All_Patients.png';
    exportgraphics(fig_traj, trajectory_filename, 'Resolution', 300);
    fprintf('Saved combined trajectory plot to: %s\n', trajectory_filename);
end
