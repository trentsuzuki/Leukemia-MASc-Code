%% Load results from a previous run

clear;
close all;

if(~isdeployed)
  cd(fileparts(matlab.desktop.editor.getActiveFilename));
end

% load("results/BacksteppingMultiStart_20241020_0331.mat");
load("results/BacksteppingMultiStart_20241021_0006.mat");

numPatients = size(patients,2);
numInitConds = size(patients(1).InitConds,2);

% Interval between plotted data points:
dataInterval = 70;

%% Time domain plots

makeTimeSeriesPlot(101,1,patients,dataInterval)
makeTimeSeriesPlot(120,1,patients,dataInterval)
makeTimeSeriesPlot(137,1,patients,dataInterval)
makeTimeSeriesPlot(176,1,patients,dataInterval)
makeTimeSeriesPlot(183,2,patients,dataInterval)

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
    
    f2 = figure('Position',[0 0 300 250],'Visible','off');
    
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

    plot(t_sim(1:dataInterval:end),x_sim(1:dataInterval:end,8),'LineWidth',1,'Color',"#0072BD")
    ylim([0 max(x_sim(:,8))])

    exportgraphics(f2,"GrantProposalPlots/SimulationResults_Patient" + patients(i).ID + ".jpg")
end