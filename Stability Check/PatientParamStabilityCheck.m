%% Stability check for patient parameter sets
% Model: Friberg (Model B from Fornari et al 2020)

clear; clc; close all;
userdir = string(java.lang.System.getProperty('user.home')) + "/";
% NOTE: The following paths are placeholders and should be updated to match your system.
codeFolder = userdir + "UofT/Research/Leukemia/Code/";
% Path to CSV file with patient parameters (Data Sheet 1 from Jost 2020 supplementary materials)
paramFilePath = codeFolder + "Jost 2020 Data Analysis/Supplementary Materials/Data Sheet 1.CSV";

%% Read patient parameters from CSV
try
    patientData = readtable(paramFilePath);
    numPatients = height(patientData);
    fprintf('Found %d patients in CSV file\n', numPatients);
catch ME
    error('Failed to read the parameter file. Please check the path: %s\nError message: %s', paramFilePath, ME.message);
end

%% Check stabilty of the Friberg model for each patient
% The stability condition is given by gamma <= gamma_star, where gamma_star
% is a function of the ratio x = k_ma / k_tr.
% We define this function first for clarity.

% Anonymous function for the stability boundary gamma_star(x)
f_gamma_star = @(x) ...
    (0.5 * (-24 - 45*x - 34*x.^2 - 9*x.^3) ./ x) + ...
    (0.5 * sqrt((576 + 2192*x + 3753*x.^2 + 3588*x.^3 + 1998*x.^4 + 612*x.^5 + 81*x.^6) ./ x.^2));

% Pre-allocate new columns in the table to store the results
patientData.gamma_star = zeros(numPatients, 1);
patientData.isStable = false(numPatients, 1);

fprintf('\n--- Stability Check Results ---\n');
fprintf('Condition for stability: Patient gamma <= calculated gamma_star\n\n');

% Loop through each patient to perform the check
for i = 1:numPatients
    % Extract necessary parameters for the current patient.
    % Note: Assumes column names 'k_tr', 'k_ma', and 'gamma' from the CSV file.
    ktr = patientData.ktr(i);
    kma = 2.375; % This value is the one used in Jost 2020
    gamma_patient = patientData.gamma(i); % The actual gamma value for the patient

    % For the drug-free steady state, the ratio Prol*/Circ* = k_ma/k_tr.
    % This is the key simplification that allows us to use the formula.
    x = kma / ktr;
    
    % Calculate the stability boundary gamma_star for this patient
    gamma_star_calculated = f_gamma_star(x);
    
    % Check if the stability condition is met (gamma <= gamma_star)
    stable = (gamma_patient <= gamma_star_calculated);
    
    % Store the results in the table
    patientData.gamma_star(i) = gamma_star_calculated;
    patientData.isStable(i) = stable;
    
    % Display the result for the current patient
    if stable
        status = 'Stable';
    else
        status = 'Unstable (prone to oscillations)';
    end
    
    % Assumes a 'Patient' column exists for identification
    fprintf('Patient %2d: gamma = %.3f, gamma_star = %.3f  -->  %s\n', ...
            patientData.ID(i), gamma_patient, gamma_star_calculated, status);
end

%% --- Final Summary ---
numStable = sum(patientData.isStable);
numUnstable = numPatients - numStable;

fprintf('\n--- Summary ---\n');
fprintf('Total patients checked: %d\n', numPatients);
fprintf('Stable parameter sets:    %d (%.1f%%)\n', numStable, (numStable/numPatients)*100);
fprintf('Unstable parameter sets:  %d (%.1f%%)\n', numUnstable, (numUnstable/numPatients)*100);

% Display a snippet of the table with the new results
fprintf('\n--- Results Table Snippet ---\n');
disp(head(patientData(:, {'ID', 'ktr', 'gamma', 'gamma_star', 'isStable'})));