% This script runs Monte Carlo simulations for the Friberg model for multiple
% patients. It reads patient-specific parameters from a CSV file and runs
% simulations for each patient with an output feedback control input
% that is updated periodically.

clear; clc; close all;
set(0,'defaultTextInterpreter','latex');

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
        error('Error: Patient data CSV not found at the specified path. Please update paramFilePath.');
    else
        rethrow(ME);
    end
end

%% Parameters
% --- Fixed Model & Control Parameters ---
k_ma = 2.375;      % Maturation rate constant (1/day)
y_r = 1.5;         % Desired steady-state neutrophil concentration (cells/L)
y_min_safe = 0.5;  % Safe lower limit in clinical practice
y_max_safe = @(Base) 4*Base; % Upper limit to be considered 'good' performance
eps = 2.1 * y_r;   % Tuning parameter for control law

% --- Simulation Parameters ---
num_simulations = 1000; % Number of simulations to run per patient
t_span = [0 200];       % Time span for each simulation (days)
sd = 0.75;              % Standard deviation for randomizing initial conditions
measurement_interval = 7;

% --- Plot options ---
neutrophil_color = [0.5 0.5 0.9 0.1];
drug_effect_color = [0.5 0.5 0.9 0.1];
y_r_color = [0.9 0.14 0.14];

%% IMPORTANT: Options to save plots as files
save_trajectory_figures = true;
save_histogram_figures = false;

%% Main Loop for Each Patient
for p_idx = 64:height(patientData) 
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
    % --- Prepare Figures for Plotting ---
    fig_traj = figure('Visible', 'off');
    set(fig_traj, 'Color', 'White', 'Theme', 'Light');
    
    % Subplot for Neutrophil Concentration
    subplot(2,1,1);
    hold on;
    title(sprintf('Patient %d: Output Feedback Control (%d-Day Updates)', patient_id, measurement_interval), ...
        'FontSize', 16);
    xlabel('Time (days)');
    ylabel('Neutrophil Conc. (y)');
    grid on;
    
    % Subplot for Drug Effect
    subplot(2,1,2);
    hold on;
    xlabel('Time (days)');
    ylabel('Drug Effect (u)');
    grid on;
    
    % Pre-allocate storage for results
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
        
        t_current = t_span(1);
        x_current = x0;
        
        T_full = []; % Stores the complete time vector
        X_full = []; % Stores the complete state trajectory
        U_full = []; % Stores the control input for plotting
        
        % Define the output feedback control law
        u_law = @(x_vec) max(0, 1 - (y_r + eps) / (x_vec(5) + eps) * (x_vec(5) / Base)^gamma);
        
        while t_current < t_span(2)
            % Calculate control input based on the state at the beginning of the interval
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
        
        % Use the aggregated results for plotting and analysis
        t = T_full;
        x = X_full;
        u = U_full;
        
        % Extract the output (neutrophil concentration)
        y = x(:, 5);
        
        % Plot the output trajectory
        subplot(2,1,1);
        if max(y) < y_max_safe(Base)
            plot(t, y, 'Color', neutrophil_color, 'HandleVisibility', 'off');
        else
            plot(t, y, 'Color', [neutrophil_color(1:3) 0.5], 'HandleVisibility', 'off');
        end
        
        % Plot the control input using stairs to show the piecewise constant nature
        subplot(2,1,2);
        t_zoh = [t(1); kron(t(2:end), [1; 1])];
        u_zoh = [kron(u(1:end-1), [1; 1]); u(end-1);];
        d = plot(t_zoh, u_zoh, '-', ...
            'Color', drug_effect_color, ...
            'Linewidth', 1, ...
            'HandleVisibility', 'off');

        if max(y) > y_max_safe(Base)
            % alpha(st, 0.5)
        end
        
        % Store results
        final_y_values(i) = y(end);
        max_y_values(i) = max(y);
        min_y_values(i) = min(y);
    end
    
    %% --- Finalize and Save Trajectory Plot ---
    subplot(2,1,1);
    yline(y_r, '--', 'Color', y_r_color, 'LineWidth', 2, 'DisplayName', 'Desired Output (y_r)');
    text(t_span(2)*0.9, y_r + 0.5, 'Reference', 'Color', y_r_color, 'FontSize', 12);
    ylim([0, max(max_y_values) * 1.1]); % Adjust y-axis
    
    subplot(2,1,2);
    ylim([0 1]); % Adjust y-axis for control plot
    
    % Save the figure and close it
    trajectory_filename = sprintf('figs/7Day_Output_fb_Patient_%d_Trajectories_Periodic.png', patient_id);
    if save_trajectory_figures
        exportgraphics(fig_traj, trajectory_filename, 'Resolution', 300);
    end
    close(fig_traj);
    
    %% Plot and Save Distribution Histograms
    fig_dist = figure('Visible', 'off');
    set(fig_dist, 'Position', [0 0 800 300])
    set(fig_dist, 'Color', 'White', 'Theme', 'Light');
    sgtitle(sprintf('Patient %d: Distribution of Simulation Outcomes (%d-Day Updates)', patient_id, measurement_interval), 'FontSize', 16);
    
    % Plot distribution of final values
    subplot(1,3,1);
    histogram(final_y_values, 30, 'HandleVisibility', 'off');
    hold on;
    xline(y_r, 'r--', 'LineWidth', 2, 'DisplayName', 'Desired (y_r)');
    title('Final Concentrations');
    xlabel('Final y');
    ylabel('Frequency');
    grid on;
    legend('show', 'interpreter', 'latex');
    
    % Plot distribution of max values
    subplot(1,3,2);
    histogram(max_y_values, 30, 'HandleVisibility', 'off');
    title('Maximum Concentrations');
    xlabel('Max y');
    grid on;
    
    % Plot distribution of min values
    subplot(1,3,3);
    histogram(min_y_values, 30, 'HandleVisibility', 'off');
    hold on;
    xline(y_min_safe, 'r--', 'LineWidth', 2, 'DisplayName', 'Safe Limit');
    title('Minimum Concentrations');
    xlabel('Min y');
    grid on;
    legend('show', 'interpreter', 'latex');
    
    % Save the figure and close it
    distribution_filename = sprintf('figs/7Day_State_fb_Patient_%d_Distributions.png', patient_id);
    if save_histogram_figures
        exportgraphics(fig_dist, distribution_filename, 'Resolution', 300);
    end
    close(fig_dist);
    
    % --- Display Statistics ---
    fprintf('  Simulation results for Patient %d:\n', patient_id);
    fprintf('  - Desired setpoint y_r: %.4f\n', y_r);
    fprintf('  - Mean of final values:   %.4f\n', mean(final_y_values));
    fprintf('  - Std Dev of final values: %.4f\n', std(final_y_values));
    
    n_safe = sum(min_y_values >= y_min_safe);
    fprintf('  - %d/%d (%.1f%%) of simulations remained above the safe limit.\n\n', ...
        n_safe, num_simulations, 100 * n_safe / num_simulations);
end

fprintf('--- All patients processed. ---\n');
