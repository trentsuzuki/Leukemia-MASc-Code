%% Load results from a previous run

clear;
close all;
% load("results/BacksteppingMultiStart_20241020_0331.mat");
load("results/BacksteppingMultiStart_20241021_0006.mat");

numPatients = size(patients,2);
numInitConds = size(patients(1).InitConds,2);

% Color map used to distinguish between different initial conditions:
colorMap = ["#0072BD","#D95319","#EDB120","#7E2F8E","#77AC30","#4DBEEE","#A2142F"];
% Interval between plotted data points:
dataInterval = 70;

%% Figure 1: Dose function

for i = 1:numPatients
    if isempty(patients(i).ID)
        continue
    end
    
    u_nom = patients(i).ControlLaw;
    
    f1 = figure('Visible','off');
    hold on
    fplot(u_nom,[0.01 30],'LineWidth',1)
    title({"Backstepping Feedback Law","Patient " + patients(i).ID})
    xlabel('ANC (G/L)')
    ylabel('6-MP Dose (mg)')
    
    exportgraphics(f1,"figures3/ControlLaw_Patient" + patients(i).ID + ".jpg")
end

%% Figure 2: Time domain plot

for i = 1:numPatients
    if isempty(patients(i).ID)
        continue
    end
    
    u_nom = patients(i).ControlLaw;

    f2 = figure('Position',[0 0 500 500],'Visible','off');
    tiledlayout(2,1,"TileSpacing","tight","Padding","tight");
    
    nexttile(1);
    hold on
    box on
    grid on
    area([0 patients(i).InitConds(1).Time(end)], [2 2], 1,"FaceColor","flat","ShowBaseLine","off","FaceAlpha",0.15,"Linestyle","None");
    title("Neutrophil Concentration for Simulated Patient " + patients(i).ID,'Interpreter','latex','FontSize',16)
    xlabel('Time (days)','Interpreter','latex')
    ylabel("ANC (G/L)",'Interpreter','latex')
    
    nexttile(2);
    hold on
    box on
    grid on
    title("6-MP Dose",'Interpreter','latex','FontSize',16)
    xlabel('Time (days)','Interpreter','latex')
    ylabel("Dose (mg)",'Interpreter','latex')

    for j = 1:numInitConds
        t_sim = patients(i).InitConds(j).Time;
        x_sim = patients(i).InitConds(j).States;

        nexttile(1);
        plot(t_sim(1:dataInterval:end),x_sim(1:dataInterval:end,8),'LineWidth',1,'Color',colorMap(j))
        
        nexttile(2);
        plot(t_sim(1:dataInterval:end),u_nom(x_sim(1:dataInterval:end,8)),'LineWidth',1,'Color',colorMap(j))
    end
    
    exportgraphics(f2,"figures3/SimulationResults_Patient" + patients(i).ID + ".jpg")
end

%% Figure 3: Phase plot

for i = 1:numPatients
    if isempty(patients(i).ID)
        continue
    end

    p = patients(i).Params;

    ubar = vpa(p.k_20*p.CL_6tgn/(p.slope*p.FM_3*p.k_me*p.F)*(1 - (patients(i).TargetANC/p.Base)^p.gamma),6);
    xbar_3 = p.FM_3*p.k_me*p.F*ubar/(p.CL_6tgn*p.k_20);
    xbar_4 = (1-p.slope*p.FM_3*p.k_me*p.F*ubar/(p.CL_6tgn*p.k_20))^(1/p.gamma)*p.Base*p.k_ma/p.k_tr;
    xbar_8 = xbar_4*p.k_tr/p.k_ma;
    xbar = [p.F*ubar/p.k_a; p.F*ubar/p.k_20; xbar_3; xbar_4; xbar_4; xbar_4; xbar_4; xbar_8];
    
    f3 = figure('Visible','off');
    hold on
    
    plot3(p.slope*xbar(3),xbar(4),xbar(8),'k.','MarkerSize',16)
    
    for j = 1:numInitConds
        x_sim = patients(i).InitConds(j).States;
        plot3(p.slope*x_sim(1:dataInterval:end,3),x_sim(1:dataInterval:end,4),x_sim(1:dataInterval:end,8),'Color',colorMap(j));
    end
    
    title({"Phase Plot","Patient " + patients(i).ID})
    xlabel('E_{drug}')
    ylabel('x_{pr}')
    zlabel('x_{ma}')
    box on
    grid minor
    view([45 45])
    
    exportgraphics(f3,"figures3/PhasePlot_Patient" + patients(i).ID + ".jpg")
end