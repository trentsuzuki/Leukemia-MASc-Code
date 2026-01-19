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

neutrophil_color = [0.5 0.5 0.9 0.1];

%% Import data

data = readtable("Supplementary Materials/Data Sheet 2.CSV");

%% Calculate Patient Statistics

measurementIndices = data.EVID == 0;
meas = data(measurementIndices,:);
dose = data(~measurementIndices,:);

patientStats = table(unique(meas.ID),'VariableNames',{'ID'});

%% Bounds on ANC

ancBounds = rowfun(@bounds,meas, ...
    "InputVariables","DV", ...
    "GroupingVariables","ID", ...
    "OutputVariableNames",["Min ANC","Max ANC"]);

patientStats = join(patientStats,ancBounds);

figure;
hold on;
plot(patientStats.("Max ANC"),".","markersize",12)
plot(patientStats.("Min ANC"),".","markersize",12)
area([0 120], [2 2], 0.5, ...
    "FaceColor","flat", ...
    "FaceAlpha",0.15, ...
    "LineStyle","none", ...
    "ShowBaseLine","off");
title("ANC Bounds for each Patient")
xlabel("Patient Number")
ylabel("ANC [G/L]")

%% Mean ANC measurement

patientStats.("Mean ANC") = rowfun(@mean,meas, ...
    "InputVariables","DV", ...
    "GroupingVariables","ID", ...
    "OutputFormat","uniform");

figure;
hold on;
plot(patientStats.("Mean ANC"),".","markersize",12)
yline(0.5,"--")
% yline(2,"--","Upper Limit of Safe Region")
area([0 120], [2 2], 0.5, ...
    "FaceColor","flat", ...
    "FaceAlpha",0.15, ...
    "LineStyle","none", ...
    "ShowBaseLine","off");
title("Mean ANC for each Patient")
xlabel("Patient Number")
ylabel("Mean ANC [G/L]")

unsafeMeanIndices = patientStats.("Mean ANC") > 2 | patientStats.("Mean ANC") < 0.5;
numPatientsWithMeanANCUnsafe = nnz(unsafeMeanIndices);
numPatients = length(unsafeMeanIndices);
percentOfPatientsWithMeanANCUnsafe = numPatientsWithMeanANCUnsafe / numPatients * 100;

%% Median ANC Measurement

patientStats.("Median ANC") = rowfun(@median,meas, ...
    "InputVariables","DV", ...
    "GroupingVariables","ID", ...
    "OutputFormat","uniform");

figure;
hold on;
plot(patientStats.("Median ANC"),".","markersize",12)
area([0 120], [2 2], 0.5, ...
    "FaceColor","flat", ...
    "FaceAlpha",0.15, ...
    "LineStyle","none", ...
    "ShowBaseLine","off");
title("Median ANC for each Patient")
xlabel("Patient Number")
ylabel("Median ANC [G/L]")

unsafeMedianIndices = patientStats.("Median ANC") > 2 | patientStats.("Median ANC") < 0.5;
numPatientsWithMedianANCUnsafe = nnz(unsafeMedianIndices);
numPatients = length(unsafeMedianIndices);
percentOfPatientsWithMedianANCUnsafe = numPatientsWithMedianANCUnsafe / numPatients * 100;

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

%% Plot all patients

figure;
hold on

for ii = 1:height(measurementGaps)
    plot(measurementGaps.daysBetweenMeasurements{ii},'o')
end

title("Days Between ANC Measurements (All Patients)")
xlabel('Measurement Index')
ylabel('Days Between Measurements')

% for ii = 1:height(measurementGaps)
%     gaps = measurementGaps.daysBetweenMeasurements{ii};
%     if any(gaps > 300)
%         disp(ii)
%     end
% end

% ANC trajectories for all patients:
figure;
hold on

for id = unique(meas.ID)'
    measForThisID = meas(meas.ID == id,:);
    plot(measForThisID.TIME, measForThisID.DV, ...
        'Color', neutrophil_color, ...
        'LineWidth', 2);
end

title("Patient ANC Trajectories")
xlabel("Time [days]")
ylabel("ANC [G/L]")

%% ANC trajectory for a single patient:
% Patient is selected by ID
% id = 118; % Large gap between measurements
id = 177; % BMI > 30, may require different dosing calculation
% id = 138;
% id = 101;
% id = 188;

measForThisID = meas(meas.ID == id,:);
doseForThisID = dose(dose.ID == id,:);

if isempty(measForThisID)
    error("No measurements found for patient ID " + id);
end

% Nominal 6MP dose calculated using patient BSA: 50 mg/m^2
nominal6MPDoseForThisID = patientStats{patientStats.ID == id,"Nominal Dose"};

% Find ANC measurements below safe lower limit
% Protocol AIEOP-BFM 2009 says dose should be lowered when this occurs
ancTooLowForThisID = measForThisID(measForThisID.DV < 0.5,:);
if isempty(ancTooLowForThisID)
    ancTooLowForThisID(1,:) = table(nan);
end

upperXlim = ceil(max(measForThisID.TIME)/50)*50;

figure;
tiledlayout(2,1);

nexttile;
hold on
plot(measForThisID.TIME, measForThisID.DV);
plot(ancTooLowForThisID.TIME, ancTooLowForThisID.DV,"*");
xline(ancTooLowForThisID.TIME,":")
area([0 upperXlim], [2 2], 0.5, ...
    "FaceColor","flat", ...
    "FaceAlpha",0.15, ...
    "LineStyle","none", ...
    "ShowBaseLine","off");

% Sanity check
% yline(patientStats{patientStats.ID == id,"Min ANC"})
% yline(patientStats{patientStats.ID == id,"Max ANC"})

title("ANC Trajectory (Patient " + id + ")")
xlabel("Time [days]")
ylabel("ANC [G/L]")
xlim([0 upperXlim])
grid on

nexttile;
hold on
stairs(doseForThisID.TIME, doseForThisID.AMT);
plot(doseForThisID.TIME, doseForThisID.AMT,'k.');
xline(ancTooLowForThisID.TIME,":")
yline(nominal6MPDoseForThisID,"--","Nominal Dose: " + round(nominal6MPDoseForThisID) + " mg")
title("6-MP Dosing (Patient " + id + ")")
xlabel("Time [days]")
ylabel("6-MP Dose [mg]")
xlim([0 upperXlim])
ylim padded
grid on
