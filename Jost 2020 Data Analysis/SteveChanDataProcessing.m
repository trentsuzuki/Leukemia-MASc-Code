%% Setup

clear;
clc;

if(~isdeployed)
  cd(fileparts(matlab.desktop.editor.getActiveFilename));
end

%% Select the desired Excel file

fileName = "PT Labs Summarized Data 20240828 (Unencrypted).xlsx";
filePath = "./data/" + fileName;

%% 

% rawData = readtable(filePath,"TextType","string");

numericVariableNames = [
    "PatientID", ...
    "Time_days_", ...
    ... "Age_years_", ...
    "Cycle", ...
    "DayInCycle", ...
    ... "Hb_g_L_", ...
    "WBC_10_9Cells_L_", ...
    ... "PLT_x10_9Cells_L_", ...
    ... "MPV_fL_", ...
    ... "RBC_x10_12Cells_L_", ...
    ... "Hct_L_L_", ...
    ... "MCV_fL_", ...
    ... "MCH_pg_", ...
    ... "MCHC_g_L_", ...
    ... "RDW_CV___", ...
    "ANC_10_9Cells_L_", ...
    ... "Lymphocytes_x10_9Cells_L_", ...
    ... "Monocytes_x10_9Cells_L_", ...
    ... "Eosinophils_x10_9Cells_L_", ...
    ... "Basophils_x10_9Cells_L_", ...
    ... "Blasts_x10_9Cells_L_", ...
    ... "nRBCs_100WBCs__100WBC_s_", ...
    ... "nRBCs__100WBC_", ...
    ... "Metamyelocytes_x10_9Cells_L_", ...
    ... "Myelocytes_x10_9Cells_L_", ...
    ... "Promyelocytes_x10_9Cells_L_", ...
    ... "ImmatureGranulocytes_x10_9Cells_L_", ...
    ... "ImmaturePlateletFraction___", ...
    ... "AtypicalLymphocytes_x10_9Cells_L_", ...
    ... "OtherCells_x10_9Cells_L_", ...
    ... "AbnormalLymphocytes_x10_9Cells_L_", ...
    ... "Prolymphocytes_x10_9Cells_L_", ...
    ... "LargeGranularLymphocytes_x10_9Cells_L_", ...
    ... "SezaryCells_x10_9Cells_L_", ...
    "x6_MPDose_mg_", ...
    "MTXDose_mg_", ...
    ... "Weight_kg_", ...
    ... "Height_m_", ...
    "BodySurfaceArea_m_2_", ...
];

nonNumericVariableNames = [
    ... "Date", ...
    "IntensificationOrContinuation", ...
    "CycleDelayedOrShortened__Y_N_", ...
    "IfY_Reason_", ...
    "ALLSubtype_BOrT_", ...
    "PhStatus_posOrNeg_"
];

selectedVariableNames = horzcat(numericVariableNames, nonNumericVariableNames)';

% First, read almost all columns as string values.
% This is done because some numbers contain extra non-standard whitespace characters which prevent proper automatic conversion to numeric datatypes.

opts = detectImportOptions(filePath);
opts = setvartype(opts, opts.VariableNames, "string");
% opts = setvartype(opts, "Date", "datetime"); % Read "Date" column as datetime

rawData = readtable(filePath, opts);

data = rawData(:,selectedVariableNames);

% Convert numeric columns to double
for varName = numericVariableNames
    data.(varName) = double(strip((data.(varName)),' '));
end

cols.ID = "PatientID";
cols.TIME = "Time_days_";
cols.AGE = "Age_years_";
cols.CYCLE = "Cycle";
cols.DAY_IN_CYCLE = "DayInCycle";
cols.BSA = "BodySurfaceArea_m_2_";
cols.TREATMENT_PHASE = "IntensificationOrContinuation";
cols.ANC = "ANC_10_9Cells_L_";
cols.WBC = "WBC_10_9Cells_L_";
cols.DOSE_6MP = "x6_MPDose_mg_";
cols.DOSE_MTX = "MTXDose_mg_";
cols.WEIGHT = "Weight_kg_";
cols.HEIGHT = "Height_m_";

%% Data filtering

% Remove rows containing all NaNs
allNanRowIndices = all(ismissing(data),2);
data = data(~allNanRowIndices,:);
disp("Removed " + nnz(allNanRowIndices) + " empty row(s)");

% Alert if there are multiple rows with the same date for each patient ID
checkAllDatesUnique = rowfun(@(x) length(unique(x)) == length(x),data, ...
    "InputVariables",cols.TIME, ...
    "GroupingVariables",cols.ID, ...
    "OutputFormat","Uniform");

disp(nnz(checkAllDatesUnique) + " patients have all unique dates in the rows of the spreadsheet")

% On remaining rows with patient ID = NaN, replace with the value in the cell above
% e.g. patient 17 has one PatientID field missing in 20240602 dataset
missingPatientIDIndices = ismissing(data.(cols.ID));
disp("Found " + nnz(missingPatientIDIndices) + " row(s) with Patient ID missing");

if nnz(missingPatientIDIndices)
    data.(cols.ID) = fillmissing(data.(cols.ID),"previous");
    disp("Filled missing Patient ID with value from previous cell");
end

%% Create struct of patient data in a struct format

clearvars patients

patientIDs = unique(data.(cols.ID));
numPatients = length(patientIDs);

patients(numPatients) = struct;

for i = 1:numPatients
    dataForThisID = data(data.(cols.ID) == patientIDs(i),:);

    patients(i).ID = patientIDs(i);
    patients(i).Time = dataForThisID.(cols.TIME);
    patients(i).ANC = dataForThisID.(cols.ANC);
    patients(i).WBC = dataForThisID.(cols.WBC);
    patients(i).Dose_6MP = dataForThisID.(cols.DOSE_6MP);
    patients(i).Dose_MTX = dataForThisID.(cols.DOSE_MTX);
    patients(i).BSA = dataForThisID.(cols.BSA);
    patients(i).Cycle = dataForThisID.(cols.CYCLE);
    patients(i).DayInCycle = dataForThisID.(cols.DAY_IN_CYCLE);
    patients(i).TreatmentPhase = dataForThisID.(cols.TREATMENT_PHASE);
end

current_time = string(datetime,'yyyyMMdd_HHmm');
save("results/ChanData_" + current_time, 'patients')

%% Create CSV file in NONMEM/nlmixr2 format

variableNames = {'ID','TIME','AMT','DV','CMT','EVID','MDV','BSA'};
nonmem_table = array2table(zeros(0,8), 'VariableNames',variableNames);

patientIDs = unique(data.(cols.ID));
DOSE_CMT = 1; % name of compartment where drug enters the model
OBS_CMT = 8; % name of the observation compartment in the model specification

for i = 1:height(data)
    dataForThisID = data(data.(cols.ID) == data{i,cols.ID},:);
    bsa = mode(dataForThisID.(cols.BSA));

    % If row contains observation, append a new row with EVID = 0
    if ~isnan(data{i,cols.ANC})
        newRow = {
            data{i,cols.ID}, ...
            data{i,cols.TIME}, ...
            0, ...
            data{i,cols.ANC}, ...
            OBS_CMT, ...
            0, ... % EVID
            0, ...
            bsa
        };
        nonmem_table(end+1,:) = newRow;
    end
    % If row contains dose, append a new row with EVID = 1
    if ~isnan(data{i,cols.DOSE_6MP})
        newRow = {
            data{i,cols.ID}, ...
            data{i,cols.TIME}, ...
            data{i,cols.DOSE_6MP}, ...
            0, ...
            DOSE_CMT, ...
            101, ... % EVID, see https://ascpt.onlinelibrary.wiley.com/doi/10.1002/psp4.12445#psp412445-sec-0006-title
            0, ...
            bsa
        };
        nonmem_table(end+1,:) = newRow;
    end
end

% Write the table to a CSV file
outputFileName = 'nonmem_data.csv';
writetable(nonmem_table, outputFileName, 'WriteVariableNames', true);

% Display confirmation
disp(['CSV file "' outputFileName '" created successfully.']);