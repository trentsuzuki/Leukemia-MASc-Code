%% Load results from a run of BacksteppingMultiStart.m

clear;
% load("results/BacksteppingMultiStart_20241020_0331.mat");
load("results/BacksteppingMultiStart_20241021_0006.mat");

numPatients = size(patients,2);
numInitConds = size(patients(1).InitConds,2);

%% Calculate max and min ANC values for each simulation run and each patient

t_start = 7; % Data from on or before this time will be ignored for min and max ANC

for i = 1:numPatients
    if isempty(patients(i).ID)
        continue
    end

    for j = 1:numInitConds
        t_sim = patients(i).InitConds(j).Time;
        x_sim = patients(i).InitConds(j).States;
        
        startInd = find(t_sim > t_start, 1);
        patients(i).InitConds(j).MinANC = min(x_sim(startInd:end,8));
        patients(i).InitConds(j).MaxANC = max(x_sim(startInd:end,8));
    end
    
    patients(i).MinANC = min([patients(i).InitConds.MinANC]);
    patients(i).MaxANC = max([patients(i).InitConds.MaxANC]);
end

%% Calculate settling time for each patient
% Settling time is determined as the largest time such that the distance
% from the target value is less than the specified tolerance

settlingTimeTol = 0.1; % Settling time tolerance, as % of target

for i = 1:numPatients
    if isempty(patients(i).ID)
        continue
    end

    for j = 1:numInitConds
        t_sim = patients(i).InitConds(j).Time;
        x_sim = patients(i).InitConds(j).States;

        notSettledInd = abs(x_sim(:,8) - patients(i).TargetANC) >= settlingTimeTol*patients(i).TargetANC;
        settleTimeInd = find(notSettledInd, 1, "last");
        patients(i).InitConds(j).SettleTime = t_sim(settleTimeInd);
    end

    patients(i).MaxSettleTime = max([patients(i).InitConds.SettleTime]);
    patients(i).MedianSettleTime = median([patients(i).InitConds.SettleTime]);
end

%% Create struct of all simulations

simulations = [patients.InitConds];

%% Plot max and min ANC statistics

maxAllowableANC = 2.5;
minAllowableANC = 0.5;

maxSafeANC = 2;
minSafeANC = 1;

h1 = figure('Position',[0 0 300 250]);
box on
hold on
histogram([simulations.MaxANC],20)
title("Max Neutrophil Concentration (G/L)","Interpreter","latex")
ylabel("Count","Interpreter","latex")
set(gca,"TickLabelInterpreter","latex")

h2 = figure('Position',[0 0 300 250]);
box on
hold on
histogram([simulations.MinANC],20)
title("Min Neutrophil Concentration (G/L)","Interpreter","latex")
ylabel("Count","Interpreter","latex")
set(gca,"TickLabelInterpreter","latex")

% exportgraphics(h1,"GrantProposalPlots/Histogram_MaxANC.jpg")
% exportgraphics(h2,"GrantProposalPlots/Histogram_MinANC.jpg")

%% Plot settling time statistics

maxSafeSettlingTime = 100;

h3 = figure('Position',[0 0 300 250]);
box on
hold on
histogram([simulations.SettleTime],20)
title("Settling Time (days)","Interpreter","latex")
ylabel("Count","Interpreter","latex")
set(gca,"TickLabelInterpreter","latex")

% exportgraphics(h3,"GrantProposalPlots/Histogram_SettlingTime.jpg")

%% Calculate fraction of patients who meet safety criteria

numSimulations = size(simulations,2);

safeSettleTime = [simulations.SettleTime] < maxSafeSettlingTime;
safeMinANC = [simulations.MinANC] > minAllowableANC;
safeMaxANC = [simulations.MaxANC] < maxAllowableANC;

fractionWithSafeSettleTime = nnz(safeSettleTime) / numSimulations;

fractionWithMinANCAboveLimit = nnz(safeMinANC) / numSimulations;
fractionWithMaxANCBelowLimit = nnz(safeMaxANC) / numSimulations;

fractionMeetingAll3Criteria = nnz(safeSettleTime & safeMinANC & safeMaxANC) / numSimulations;
