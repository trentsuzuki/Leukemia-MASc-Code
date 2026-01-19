%% Purpose
% This script processes patient data from the supplementary materials
% provided in Jost et al., 2020, "Model-Based Simulation of Maintenance
% Therapy of Childhood Acute Lymphoblastic Leukemia"

%% Setup

clear;
clc;

if(~isdeployed)
  cd(fileparts(matlab.desktop.editor.getActiveFilename));
end

neutrophil_color = [0.5 0.5 0.9 1];
drug_effect_color = [0.5 0.5 0.9 1];

%% Import data

data = readtable("Supplementary Materials/Data Sheet 2.CSV");

%% Calculate Patient Statistics

measurementIndices = data.EVID == 0;
meas = data(measurementIndices,:);
dose = data(~measurementIndices,:);

patientStats = table(unique(meas.ID),'VariableNames',{'ID'});

%% Check patient body measurements are constant over the data set

uniqueWeightCheck = rowfun(@(weight) length(unique(weight)) <= 1,meas, ...
    "InputVariables","WEIGHT", ...
    "GroupingVariables","ID", ...
    "OutputVariableNames","isWeightUnique");

uniqueHeightCheck = rowfun(@(height) length(unique(height)) <= 1,meas, ...
    "InputVariables","HEIGHT", ...
    "GroupingVariables","ID", ...
    "OutputVariableNames","isHeightUnique");

uniqueBSACheck = rowfun(@(bsa) length(unique(bsa)) <= 1,meas, ...
    "InputVariables","BSA", ...
    "GroupingVariables","ID", ...
    "OutputVariableNames","isBSAUnique");

if any(uniqueWeightCheck.isWeightUnique == 0) ...
        || any(uniqueHeightCheck.isHeightUnique == 0) ...
        || any(uniqueBSACheck.isBSAUnique == 0)
    error("Patient body measurements changed when they were expected to be constant")
end

%% Calculate patient statistics

patientStats.("WEIGHT") = rowfun(@(arr) arr(1),meas, ...
    "InputVariables","WEIGHT", ...
    "GroupingVariables","ID", ...
    "OutputFormat","uniform");
patientStats.("HEIGHT") = rowfun(@(arr) arr(1),meas, ...
    "InputVariables","HEIGHT", ...
    "GroupingVariables","ID", ...
    "OutputFormat","uniform");
patientStats.("BSA") = rowfun(@(arr) arr(1),meas, ...
    "InputVariables","BSA", ...
    "GroupingVariables","ID", ...
    "OutputFormat","uniform");

% Calculate patient BMI:
% The protocol mentions the possibility for dose adjustment for obese patients (BMI >= 30)
patientStats.("BMI") = patientStats.WEIGHT ./ (patientStats.HEIGHT/100).^2;

obesePatientIDs = patientStats.ID(patientStats.BMI >= 30,:);

% Calculate nominal dose for each patient:
patientStats.("Nominal Dose") = 50*patientStats.BSA;

% Time between measurements:
measurementGaps = rowfun(@(measTime) {diff(measTime)},meas, ...
    "InputVariables","TIME", ...
    "GroupingVariables","ID", ...
    "OutputVariableNames","daysBetweenMeasurements");

%% Plot all patients in a loop

for id = unique(meas.ID)'
    fig_traj = figure('Visible', 'off');
    set(fig_traj, 'Color', 'White', 'Theme', 'Light');
    measForThisID = meas(meas.ID == id,:);
    doseForThisID = dose(dose.ID == id,:);

    % Subplot for Neutrophil Concentration
    subplot(2,1,1);
    hold on;
    xlabel('Time (days)');
    ylabel('ANC (G/L)');
    grid on;
    plot(measForThisID.TIME, measForThisID.DV, ...
        '*-', ...
        'Color', neutrophil_color, ...
        'LineWidth', 2);
    title(['Patient ', num2str(id), ' Treatment Outcome'], ...
    'FontSize', 16);
    
    % Subplot for Drug Effect
    subplot(2,1,2);
    hold on;
    xlabel('Time (days)');
    ylabel('Daily 6-MP Dose (mg)');
    grid on;
    stairs(doseForThisID.TIME, doseForThisID.AMT, ...
        'Color', drug_effect_color, ...
        'Linewidth', 2);
    ylim([0 150])

    trajectory_filename = sprintf('figs/Jost2020_Patient_%d_Outcome.png', id);
    exportgraphics(fig_traj, trajectory_filename, 'Resolution', 300);
    
    close(fig_traj)
end
