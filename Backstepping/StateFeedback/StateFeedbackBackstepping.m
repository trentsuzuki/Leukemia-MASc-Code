%% Setup

clear;clc;

cd(fileparts(matlab.desktop.editor.getActiveFilename));
topLevelDirectory = "C:\Users\trent\UofT\Research\Leukemia\Code";
addpath(topLevelDirectory);

%% Create object holding symbolic expressions for the model dynamics

syms x(t) [8 1]
syms F k_a k_20 FM_3 k_me CL_6tgn k_ma
syms Base k_tr gamma slope
syms u

% First column: x dot
dxdt(1,1) = diff(x1,t);
dxdt(2,1) = diff(x2,t);
dxdt(3,1) = diff(x3,t);
dxdt(4,1) = diff(x4,t);
dxdt(5,1) = diff(x5,t);
dxdt(6,1) = diff(x6,t);
dxdt(7,1) = diff(x7,t);
dxdt(8,1) = diff(x8,t);
% Second column: Dynamics for x dot
dxdt(1,2) = -k_a*x1 + F*u;
dxdt(2,2) = k_a*x1 - k_20*x2;
dxdt(3,2) = FM_3*k_me*x2 - CL_6tgn*x3;
dxdt(4,2) = k_tr*x4*(1-slope*x3)*(Base/x8)^gamma - k_tr*x4;
dxdt(5,2) = k_tr*(x4 - x5);
dxdt(6,2) = k_tr*(x5 - x6);
dxdt(7,2) = k_tr*(x6 - x7);
dxdt(8,2) = k_tr*x7 - k_ma*x8;

%% Write expressions for controller using backstepping

syms x_d
syms C [8 1]
syms a [8 1]

z8 = x8 - x_d;

alpha7 = 1/k_tr*((k_ma-C1)*z8 + k_ma*x_d);
z7 = x7 - alpha7;

% Easiest approach: alpha4,5,6 = alpha7
% This approach ensures positivity of control laws alpha_i for i=4,5,6,7
alpha6 = alpha7;
alpha5 = alpha6;
alpha4 = alpha5;
z6 = x6 - alpha6;
z5 = x5 - alpha5;
z4 = x4 - alpha4;

% Complete backstepping: Individually-tuned control laws for all states 
% alpha6 = (1-C2/k_tr)*z7 + alpha7 + 1/k_tr*diff(alpha7,t) - z8;
% alpha6 = subs(alpha6, dxdt(:,1), dxdt(:,2));
% z6 = x6 - alpha6;
% 
% alpha5 = (1-C3/k_tr)*z6 + alpha6 + 1/k_tr*diff(alpha6,t) - z7;
% alpha5 = subs(alpha5, dxdt(:,1), dxdt(:,2));
% z5 = x5 - alpha5;
% 
% alpha4 = (1-C4/k_tr)*z5 + alpha5 + 1/k_tr*diff(alpha5,t) - z6;
% alpha4 = subs(alpha4, dxdt(:,1), dxdt(:,2));
% z4 = x4 - alpha4;

alpha3 = 1/slope*(1 - (x8/Base)^gamma * (1 + (diff(alpha4,t) - C5*z4)/(k_tr*x4) - z5/x4));
alpha3 = subs(alpha3, dxdt(:,1), dxdt(:,2));
z3 = x3 - alpha3;

alpha2 = 1/(FM_3*k_me)*((CL_6tgn-C6)*z3 + CL_6tgn*alpha3 + diff(alpha3,t) + k_tr*slope*z4*x4*(Base/x8)^gamma);
alpha2 = subs(alpha2, dxdt(:,1), dxdt(:,2));
z2 = x2 - alpha2;

alpha1 = 1/k_a*((k_20-C7)*z2 + k_20*alpha2 + diff(alpha2,t) - FM_3*k_me*z3);
alpha1 = subs(alpha1, dxdt(:,1), dxdt(:,2));
z1 = x1 - alpha1;

u = 1/F*((k_a-C8)*z1 + k_a*alpha1 + diff(alpha1,t) - k_a*z2);
u = subs(u, dxdt(:,1), dxdt(:,2));

%% Set patient parameters

patientID = 101;
BSA = 1.9833;
p.Base = 2.5847;
p.k_tr = 0.15455;
p.gamma = 0.975662;
p.slope = 0.24489;

patientID = 120;
BSA = 0.6491;
p.Base = 2.7738;
p.k_tr = 0.078438;
p.gamma = 1.2359;
p.slope = 0.59982;

patientID = 165;
BSA = 0.61547;
p.Base = 2.711;
p.k_tr = 0.14186;
p.gamma = 0.80269;
p.slope = 0.094133;
 
patientID = 183;
BSA = 0.58481;
p.Base = 2.5808;
p.k_tr = 0.14815;
p.gamma = 1.1076;
p.slope = 0.24548;

p.F = 0.22;
p.k_a = 31.2;
p.k_20 = 12.72;
p.FM_3 = 0.019;
p.k_me = 9.9216;
p.CL_6tgn = 0.219*BSA^1.16; % From Table 2. Not Equation (2), which uses incorrect units.
p.k_ma = 2.3765;

symbols(1,:) = [F, p.F];
symbols(2,:) = [k_a, p.k_a];
symbols(3,:) = [k_20, p.k_20];
symbols(4,:) = [FM_3, p.FM_3];
symbols(5,:) = [k_me, p.k_me];
symbols(6,:) = [CL_6tgn, p.CL_6tgn];
symbols(7,:) = [k_ma, p.k_ma];
% Patient-specific parameters:
symbols(8,:) = [Base, p.Base];
symbols(9,:) = [k_tr, p.k_tr];
symbols(10,:) = [slope, p.slope];
symbols(11,:) = [gamma, p.gamma];
symbols(12,:) = [x_d, 1.5];

% 0 < C1 <= k_ma
C_values(:,1) = [C1;C2;C3;C4;C5;C6;C7;C8];
C_values(:,2) = [p.k_ma/2;.01;.01;.1;10;10;10;10];

u_nom = subs(u, symbols(:,1), symbols(:,2));
u_nom = subs(u_nom,C_values(:,1),C_values(:,2));
u_nom(a) = subs(u_nom,x,a);
u_fun = matlabFunction(u_nom,'Vars',{a});

% Constrain the control law to the range [0,250]
u_sat = @(x) max(min(u_fun(x),250),0);

% figure;
% fplot(u_nom,[0.01 5])
% x_values = 0.05:0.05:5;
% % plot(x_values,u_nom(x_values))
% title('Backstepping Feedback Law')
% xlabel('ANC (G/L)')
% ylabel('6-MP Dose (mg)')

x0 = [0; 0; 0; (p.Base*p.k_ma)/p.k_tr; (p.Base*p.k_ma)/p.k_tr; (p.Base*p.k_ma)/p.k_tr; (p.Base*p.k_ma)/p.k_tr; p.Base];

% Evaluate control law at initial condition:
disp("u(x(0)): " + double(u_sat(x0)) + " mg")

%% Simulation

[t_sim,x_sim] = ode45(@(t,x) JostModelWithParams(t,x,u_sat(x),p), 0:1:600, x0);

%% Lyapunov function:

% Quadratic Lyapunov function:
V = 1/2*(C8*z1^2 + C7*z2^2 + C6*z3^2 + C5*z4^2 + C4*z5^2 + C3*z6^2 + C2*z7^2 + C1*z8^2);
V = subs(V,symbols(:,1),symbols(:,2));
V = subs(V,C_values(:,1),C_values(:,2));
V(a) = subs(V,x,a);
V_fun = matlabFunction(V,'Vars',{a});

% Total time derivative of Lyapunov function:
% V_dot = subs(diff(V,t),dxdt(:,1),dxdt(:,2));
% V_dot = subs(V_dot,symbols(:,1),symbols(:,2));
% V_dot = subs(V_dot,C_values(:,1),C_values(:,2));
% V_dot(a) = subs(V_dot,x,a);
% V_dot_fun = matlabFunction(V_dot,'Vars',{a});
% vpa(V_dot(x0))

figure;
plot(t_sim,log(V_fun(x_sim')),'linewidth',1,'Color',"#0072BD")
title("Lyapunov Function")
xlabel("Time (days)")
ylabel("V(x)")

%% Time domain plot:

f = figure('Position',[0 0 500 500]);
tiledlayout(2,1,"TileSpacing","tight","Padding","tight");

nexttile;
hold on
box on
grid on
area([0 t_sim(end)], [2 2], 0.5,"FaceColor","flat","ShowBaseLine","off","FaceAlpha",0.15,"Linestyle","None");
plot(t_sim,x_sim(:,8),'linewidth',1,'Color',"#0072BD")
% yline(x_target,'k--')
title({"Neutrophil Concentration", "for Simulated Patient " + patientID},'Interpreter','latex','FontSize',16)
xlabel('Time (days)','Interpreter','latex')
ylabel("ANC (G/L)",'Interpreter','latex')

nexttile;
hold on
box on
grid on
plot(t_sim,u_sat(x_sim'),'linewidth',1)
title("6-MP Dose",'Interpreter','latex','FontSize',16)
xlabel('Time (days)','Interpreter','latex')
ylabel("Dose (mg)",'Interpreter','latex')

exportgraphics(f,"figures/SaturatedSimulationResults_Patient" + patientID + ".jpg")

%% Phase plot:

ubar = double(vpa(p.k_20*p.CL_6tgn/(p.slope*p.FM_3*p.k_me*p.F)*(1 - (symbols(12,2)/p.Base)^p.gamma),6));
disp("Equilibrium drug dose: " + ubar + " mg")

% Equilibrium point:
xbar_4 = (1-p.slope*p.FM_3*p.k_me*p.F*ubar/(p.CL_6tgn*p.k_20))^(1/p.gamma)*(p.Base*p.k_ma/p.k_tr);
xbar = [p.F*ubar/p.k_a; 
    p.F*ubar/p.k_20; 
    p.FM_3*p.k_me*p.F*ubar/(p.CL_6tgn*p.k_20); 
    xbar_4; 
    xbar_4; 
    xbar_4; 
    xbar_4; 
    (1-p.slope*p.FM_3*p.k_me*p.F*ubar/(p.CL_6tgn*p.k_20))^(1/p.gamma)*p.Base];

figure;
hold on
plot3(p.slope*x_sim(:,3),x_sim(:,4),x_sim(:,8));

E_vals = p.slope*min(x_sim(:,3)):0.01:p.slope*max(x_sim(:,3))+0.02;
% plot3(E_vals,repmat(xbar(4),size(E_vals)),p.Base*(1-E_vals).^(1/p.gamma))
plot3(p.slope*xbar(3),xbar(4),xbar(8),'k.','MarkerSize',16)

points(1,:) = [E_vals(1) x_sim(1,4) x_sim(1,8)];
points(2,:) = [E_vals(1) 0 x_sim(1,8)];
points(3,:) = [E_vals(end) 0 p.Base*(1-E_vals(end))^(1/p.gamma)];
points(4,:) = [E_vals(end) x_sim(1,4) p.Base*(1-E_vals(end))^(1/p.gamma)];
% fill3(points(:,1),points(:,2),points(:,3),'r','FaceAlpha',0.3)

title("Phase Plot")
xlabel('E_{drug}')
ylabel('x_{pr}')
zlabel('x_{ma}')
box on
grid minor
view([45 45])

