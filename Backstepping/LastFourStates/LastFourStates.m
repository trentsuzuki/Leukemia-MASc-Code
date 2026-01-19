%% Setup

clear;
clc;
cd(fileparts(matlab.desktop.editor.getActiveFilename));

topLevelDir = "C:\Users\trent\UofT\Research\Leukemia\Code";
addpath(topLevelDir)

%% Symbolic functions

syms x8 x_d C1
syms F k_a k_20 FM_3 k_me CL_6tgn k_ma
syms Base k_tr gamma slope

alpha_4(x8) = (k_ma-C1)/k_tr*x8 + C1/k_tr*x_d;

%% Import patient-specific parameter values from Jost 2020 Supplementary material

% Data Sheet 1
opts1 = delimitedTextImportOptions("NumVariables", 5);

% Specify range and delimiter
opts1.DataLines = [2, Inf];
opts1.Delimiter = ",";

% Specify column names and types
opts1.VariableNames = ["ID", "Base", "ktr", "slope", "gamma"];
opts1.VariableTypes = ["double", "double", "double", "double", "double"];

% Specify file level properties
opts1.ExtraColumnsRule = "ignore";
opts1.EmptyLineRule = "read";

dataSheet1Path = topLevelDir + "\Jost 2020 Data Analysis\Supplementary Materials\Data Sheet 1.CSV";

patientSpecificParams = readtable(dataSheet1Path, opts1);

% Data Sheet 2
opts2 = delimitedTextImportOptions("NumVariables", 11);

% Specify range and delimiter
opts2.DataLines = [1, Inf];
opts2.Delimiter = ",";

% Specify column names and types
opts2.VariableNames = ["ID", "TIME", "DV", "CMT", "AMT", "MDV", "EVID", "IPRED", "WEIGHT", "HEIGHT", "BSA"];
opts2.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];

% Specify file level properties
opts2.ExtraColumnsRule = "ignore";
opts2.EmptyLineRule = "read";

dataSheet2Path = topLevelDir + "\Jost 2020 Data Analysis\Supplementary Materials\Data Sheet 2.CSV";

DataSheet2 = readtable(dataSheet2Path, opts2);

%% Initial Conditions

numInitConds = 5;
randomStateUpperLims = [60;60;60;7];
% initConds = randomStateUpperLims.*rand(4,numInitConds);

initConds = [
    7.8025	5.7521	6.4912	6.8678	4.8679;
    3.8974	0.5978	7.3172	1.8351	4.3586;
    2.4169	2.3478	6.4775	3.6848	4.4678;
    24.2347	21.1895	27.0554	37.5371	18.3810;
    5.7873	49.2716	32.8205	46.8136	30.5105;
    7.9184	0.9242	17.7792	4.8675	30.6463;
    56.5230	2.5814	44.6816	55.7632	49.0577;
    6.6929	1.1829	1.3227	5.4300	5.5638;
];

% Replace first four states with NaN
initConds(1:4,:) = nan(4,numInitConds);

disp(initConds);

%% Simulate patients from Jost 2020

clearvars patients

numPatients = height(patientSpecificParams);
patients(numPatients) = struct;

for i = 1:numPatients
    innerLoopStartTime = datetime("now");

    patients(i).ID = patientSpecificParams.ID(i);
    patients(i).InitConds(numInitConds) = struct;
    
    % Patient-specific parameters
    params.Base = patientSpecificParams.Base(i);
    params.k_tr = patientSpecificParams.ktr(i);
    params.gamma = patientSpecificParams.gamma(i);
    params.slope = patientSpecificParams.slope(i);
    params.BSA = mode(DataSheet2.BSA(DataSheet2.ID == patients(i).ID));
    
    % Population Parameters
    params.F = 0.22;
    params.k_a = 31.2;
    params.k_20 = 12.72;
    params.FM_3 = 0.019;
    params.k_me = 9.9216;
    params.CL_6tgn = 0.219*params.BSA^1.16; % From Table 2. Not Equation (2), which uses incorrect units.
    params.k_ma = 2.3765;

    % Controller parameters
    params.C1 = params.k_ma; %0.08;
    params.x_d = 1.5;

    patients(i).Params = params;

    alpha_4_num(x8) = vpa(subs(max(alpha_4(x8),0), fieldnames(params), struct2cell(params)),3);

    patients(i).TargetANC = params.x_d;
    patients(i).ControlLaw = alpha_4_num;
    
    % Simulation
    Tmax = 200;
    
    for j = 1:numInitConds
        x0 = initConds(:,j);
        
        [t_sim,x_sim] = ode45(@(t,x) JostModelLastFourStates(t,x,alpha_4_num(x(8)),params), [0 Tmax], x0);
        
        patients(i).InitConds(j).Time = t_sim;
        patients(i).InitConds(j).States = x_sim;
    end

    innerLoopStopTime = datetime("now");
    patients(i).RunTime = innerLoopStopTime - innerLoopStartTime;

    disp("Completed simulations for patient " + patients(i).ID + " (" + i + "/" + numPatients + ")");
end

current_time = string(datetime,'yyyyMMdd_HHmm');

save("../results/LastFourStates_" + current_time, 'patients')
