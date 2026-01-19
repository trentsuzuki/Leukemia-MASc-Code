%% Load results from a previous run

clear;
close all;

cd(fileparts(matlab.desktop.editor.getActiveFilename));
load("../results/LastFourStates_20241202_1925.mat");

numPatients = size(patients,2);
numInitConds = size(patients(1).InitConds,2);

dataPlotInterval = 1;

%% Check final neutrophil concentrations are near target

convergenceCheck = [];
tol = 0.01;

for i = 1:size(patients,2)
    for j = 1:size(patients(i).InitConds,1)
        convergenceCheck(end+1) = patients(i).InitConds(j).States(end,8) - patients(i).Params.x_d < tol;
    end
end

disp("All simulations converged to target concentration?: " + all(convergenceCheck))

%% Time domain plots

for patientID = [101,120,137,176,183]
    for initCond = [1,2,3,4,5]
        makeTimeSeriesPlot(patientID,initCond,patients,dataPlotInterval)
    end
end

%% Helper functions

function makeTimeSeriesPlot(id,initCond,patients,dataInterval)
    i = find([patients.ID] == id);
    j = initCond;

    numPatients = size(patients,2);
    numInitConds = size(patients(1).InitConds,2);

    if i < 1
        return;
    end
    
    u_nom = patients(i).ControlLaw;
    
    f2 = figure('Position',[0 0 300 250]);
    
    hold on
    box on
    grid on
    area([0 patients(i).InitConds(1).Time(end)], [2 2], 1,"FaceColor","flat","ShowBaseLine","off","FaceAlpha",0.15,"Linestyle","None");
    title({"Simulated Patient " + patients(i).ID},'Interpreter','latex')
    xlabel('Time (days)','Interpreter','latex')
    ylabel("Neutrophil Concentration (G/L)",'Interpreter','latex')
    set(gca,'TickLabelInterpreter','latex')
    
    t_sim = patients(i).InitConds(j).Time;
    x_sim = patients(i).InitConds(j).States;

    plot(t_sim(1:dataInterval:end),x_sim(1:dataInterval:end,:),'LineWidth',1,'Color',"#0072BD")
    ylim([0 max(x_sim(:,8))])

    % exportgraphics(f2,"GrantProposalPlots/SimulationResults_Patient" + patients(i).ID + ".jpg")
end