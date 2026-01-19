%% Robust control of the model from Jost 2020
% State feedback using H2 synthesis

%% Setup
clear;
clc;
close all;

if(~isdeployed)
  cd(fileparts(matlab.desktop.editor.getActiveFilename));
end

%% Parameters

n = 8; % number of states

% Population parameters
BSA = 0.82;
F = 0.22;
k_a = 31.2;
k_20 = 12.72;
FM_3 = 0.019;
k_me = 9.9216;
CL_6tgn = 0.219*BSA^1.16;
k_ma = 2.3765;

% Bounds for patient specific parameters

% Ideal bounds to cover all patients from Jost 2020:
% pmin.Base = 1.5;
% pmax.Base = 5;
% pmin.k_tr = 0.05;
% pmax.k_tr = 0.25;
% pmin.gamma = 0.5;
% pmax.gamma = 1.3;
% pmax.slope = 0.7;
% pmin.slope = 0.08;

% Bounds over which LMIs are feasible:
pmin.Base = 2.2;
pmax.Base = 2.4;
pmin.k_tr = 0.15;
pmax.k_tr = 0.16;
pmin.gamma = 0.75;
pmax.gamma = 0.8;
pmin.slope = 0.22;
pmax.slope = 0.3;

% Patient-specific parameters for sample patient (Table 5)
Base_sample = 2.34;
k_tr_sample = 0.148;
slope_sample = 0.242;
gamma_sample = 0.769;

% Initial condition for sample patient
x0_jost = [0; 0; 0; (Base_sample*k_ma)/k_tr_sample; (Base_sample*k_ma)/k_tr_sample; (Base_sample*k_ma)/k_tr_sample; (Base_sample*k_ma)/k_tr_sample; Base_sample];

%% Equilibrium

xbar_8 = 1.5; % Desired ANC

% Calculate average equilibrium dose (ubar) over the uncertainty set
expr1 = @(b) ((xbar_8./b).^pmax.gamma - (xbar_8./b).^pmin.gamma) .* log(xbar_8 ./ b);
ubar_integral = integral(expr1, pmin.Base, pmax.Base);
ubar_ave = CL_6tgn*k_20/(FM_3*k_me*F*(pmax.slope-pmin.slope)) * log(pmax.slope/pmin.slope) * (1 - 1/(pmax.Base-pmin.Base)/(pmax.gamma-pmin.gamma) * ubar_integral);

xbar_4 = k_ma*xbar_8/(pmax.k_tr-pmin.k_tr) * log(pmax.k_tr/pmin.k_tr);
xbar_3 = FM_3*k_me*F*ubar_ave / (CL_6tgn*k_20);

xbar_ave = [F*ubar_ave/k_a; F*ubar_ave/k_20; xbar_3; xbar_4; xbar_4; xbar_4; xbar_4; xbar_8];

%% Representation of Uncertain System
% A(delta) is affine in delta
% delta belongs to a box-like set

A0 = zeros(n);
A0(1,1) = -k_a;
A0(2,1) = k_a;
A0(2,2) = -k_20;
A0(3,2) = FM_3*k_me;
A0(3,3) = -CL_6tgn;
A0(8,8) = -k_ma;

A1 = zeros(n);
A1(4,4) = -1;
A1(5,4) = 1;
A1(5,5) = -1;
A1(6,5) = 1;
A1(6,6) = -1;
A1(7,6) = 1;
A1(7,7) = -1;
A1(8,7) = 1;

A2 = zeros(n);
A2(4,4) = 1;

A3 = zeros(n);
A3(4,8) = -xbar_4/xbar_8;

A4 = zeros(n);
A4(4,8) = xbar_4*xbar_3/xbar_8;

A5 = zeros(n);
A5(4,3) = -xbar_4;

% Affine (w.r.t. uncertainty) representation of dynamics
A = @(delta) A0 + A1*delta(1) + A2*delta(2) + A3*delta(3) + A4*delta(4) + A5*delta(5);

delta_min = [pmin.k_tr, pmin.k_tr*(pmin.Base/xbar_8)^pmin.gamma, pmin.k_tr*pmin.gamma*(pmin.Base/xbar_8)^pmin.gamma, pmin.k_tr*pmin.gamma*(pmin.Base/xbar_8)^pmin.gamma*pmin.slope, pmin.k_tr*(pmin.Base/xbar_8)^pmin.gamma*pmin.slope];
delta_max = [pmax.k_tr, pmax.k_tr*(pmax.Base/xbar_8)^pmax.gamma, pmax.k_tr*pmax.gamma*(pmax.Base/xbar_8)^pmax.gamma, pmax.k_tr*pmax.gamma*(pmax.Base/xbar_8)^pmax.gamma*pmax.slope, pmax.k_tr*(pmax.Base/xbar_8)^pmax.gamma*pmax.slope];

theta_to_delta = @(theta) [theta(2), theta(2)*(theta(1)/xbar_8)^theta(4), theta(2)*theta(4)*(theta(1)/xbar_8)^theta(4), theta(2)*theta(4)*(theta(1)/xbar_8)^theta(4)*theta(3), theta(2)*(theta(1)/xbar_8)^theta(4)*theta(3)];
delta_sample = theta_to_delta([Base_sample;k_tr_sample;slope_sample;gamma_sample]);
A_nom = A(delta_sample);

num_dimensions = length(delta_min);
num_corners = 2^num_dimensions;
corners = zeros(num_corners, num_dimensions);
A_corners = cell(num_corners,1);

% Generate all corners of the convex set
for i = 0:num_corners-1
    binary_index = dec2bin(i, num_dimensions) - '0'; % Convert to binary array
    corners(i+1, :) = delta_min .* (1 - binary_index) + delta_max .* binary_index;
    A_corners{i+1} = A(corners(i+1,:));
    eig(A_corners{i+1})
end

Bu = [F; zeros(7,1)];
Cy = [zeros(1,7) 1];

sys_nom = compreal(ss(A_nom,Bu,Cy,[]),"o"); % observable canonical form

%% Robust H2 synthesis (LQR with initial condition as impulsive disturbance)

m = 1;
nz = 2;
nw = 1;

Bw = x0_jost - xbar_ave;
Cz = [zeros(1,n-1) 1e-3; zeros(m,8)];
Dzw = zeros(nz,nw);
Dzu = [zeros(1,m); 1e-3];

X = sdpvar(n);
Z = sdpvar(m,n,'full');
W = sdpvar(nz);
small = 1e-5;

Constraints = [
    X >= small*eye(n);
    [X,(Cz*X+Dzu*Z)';(Cz*X+Dzu*Z),W] >= small*eye(n+nz);
    [A_corners{1},Bu]*[X;Z] + ([A_corners{1},Bu]*[X;Z])' + Bw*Bw' <= -small*eye(n);
    [A_corners{2},Bu]*[X;Z] + ([A_corners{2},Bu]*[X;Z])' + Bw*Bw' <= -small*eye(n);
    [A_corners{3},Bu]*[X;Z] + ([A_corners{3},Bu]*[X;Z])' + Bw*Bw' <= -small*eye(n);
    [A_corners{4},Bu]*[X;Z] + ([A_corners{4},Bu]*[X;Z])' + Bw*Bw' <= -small*eye(n);
    [A_corners{5},Bu]*[X;Z] + ([A_corners{5},Bu]*[X;Z])' + Bw*Bw' <= -small*eye(n);
    [A_corners{6},Bu]*[X;Z] + ([A_corners{6},Bu]*[X;Z])' + Bw*Bw' <= -small*eye(n);
    [A_corners{7},Bu]*[X;Z] + ([A_corners{7},Bu]*[X;Z])' + Bw*Bw' <= -small*eye(n);
    [A_corners{8},Bu]*[X;Z] + ([A_corners{8},Bu]*[X;Z])' + Bw*Bw' <= -small*eye(n);
    [A_corners{9},Bu]*[X;Z] + ([A_corners{9},Bu]*[X;Z])' + Bw*Bw' <= -small*eye(n);
    [A_corners{10},Bu]*[X;Z] + ([A_corners{10},Bu]*[X;Z])' + Bw*Bw' <= -small*eye(n);
    [A_corners{11},Bu]*[X;Z] + ([A_corners{11},Bu]*[X;Z])' + Bw*Bw' <= -small*eye(n);
    [A_corners{12},Bu]*[X;Z] + ([A_corners{12},Bu]*[X;Z])' + Bw*Bw' <= -small*eye(n);
    [A_corners{13},Bu]*[X;Z] + ([A_corners{13},Bu]*[X;Z])' + Bw*Bw' <= -small*eye(n);
    [A_corners{14},Bu]*[X;Z] + ([A_corners{14},Bu]*[X;Z])' + Bw*Bw' <= -small*eye(n);
    [A_corners{15},Bu]*[X;Z] + ([A_corners{15},Bu]*[X;Z])' + Bw*Bw' <= -small*eye(n);
    [A_corners{16},Bu]*[X;Z] + ([A_corners{16},Bu]*[X;Z])' + Bw*Bw' <= -small*eye(n);
    [A_corners{17},Bu]*[X;Z] + ([A_corners{17},Bu]*[X;Z])' + Bw*Bw' <= -small*eye(n);
    [A_corners{18},Bu]*[X;Z] + ([A_corners{18},Bu]*[X;Z])' + Bw*Bw' <= -small*eye(n);
    [A_corners{19},Bu]*[X;Z] + ([A_corners{19},Bu]*[X;Z])' + Bw*Bw' <= -small*eye(n);
    [A_corners{20},Bu]*[X;Z] + ([A_corners{20},Bu]*[X;Z])' + Bw*Bw' <= -small*eye(n);
    [A_corners{21},Bu]*[X;Z] + ([A_corners{21},Bu]*[X;Z])' + Bw*Bw' <= -small*eye(n);
    [A_corners{22},Bu]*[X;Z] + ([A_corners{22},Bu]*[X;Z])' + Bw*Bw' <= -small*eye(n);
    [A_corners{23},Bu]*[X;Z] + ([A_corners{23},Bu]*[X;Z])' + Bw*Bw' <= -small*eye(n);
    [A_corners{24},Bu]*[X;Z] + ([A_corners{24},Bu]*[X;Z])' + Bw*Bw' <= -small*eye(n);
    [A_corners{25},Bu]*[X;Z] + ([A_corners{25},Bu]*[X;Z])' + Bw*Bw' <= -small*eye(n);
    [A_corners{26},Bu]*[X;Z] + ([A_corners{26},Bu]*[X;Z])' + Bw*Bw' <= -small*eye(n);
    [A_corners{27},Bu]*[X;Z] + ([A_corners{27},Bu]*[X;Z])' + Bw*Bw' <= -small*eye(n);
    [A_corners{28},Bu]*[X;Z] + ([A_corners{28},Bu]*[X;Z])' + Bw*Bw' <= -small*eye(n);
    [A_corners{29},Bu]*[X;Z] + ([A_corners{29},Bu]*[X;Z])' + Bw*Bw' <= -small*eye(n);
    [A_corners{30},Bu]*[X;Z] + ([A_corners{30},Bu]*[X;Z])' + Bw*Bw' <= -small*eye(n);
    [A_corners{31},Bu]*[X;Z] + ([A_corners{31},Bu]*[X;Z])' + Bw*Bw' <= -small*eye(n);
    [A_corners{32},Bu]*[X;Z] + ([A_corners{32},Bu]*[X;Z])' + Bw*Bw' <= -small*eye(n);
];

Cost = trace(W);
options = sdpsettings('solver','sdpt3','verbose',1,'sdpt3.maxit',1000);
sol = optimize(Constraints,Cost,options);

sol.info
K_H2 = value(Z)*inv(value(X));
sqrt(value(Cost))

%% Simulate H2 controller on linear system

t = 0:0.1:400;
w = zeros(length(t),nw);

x_H2 = cell(num_corners,1);
u_H2 = nan(num_corners,length(t));
output_names = {'x_1','x_2','x_3','x_4','x_5','x_6','x_7','x_8'};

for i = 1:num_corners
    sys_i = ss(A(corners(i,:))+Bu*K_H2,zeros(n,nw),eye(n),zeros(n,nw),'TimeUnit','Days','OutputName',output_names);
    x_H2(i) = {lsim(sys_i,w,t,x0_jost-xbar_ave)};
    u_H2(i,:) = K_H2*x_H2{i}' + ubar_ave;
end

f_H2 = figure('Position',[0 0 800 500]);
tiledlayout(2,1,"TileSpacing","tight","Padding","compact");

nexttile;
hold on

for i = 1:num_corners
    x = x_H2{i};
    plot(t,x(:,8) + xbar_8,'linewidth',2)
end

yline(xbar_8,'k--')
title("Disturbance Impulse Response (Neutrophil Count in G/L)",'Interpreter','latex','FontSize',16)
xlabel('Time (days)','Interpreter','latex')
set(gca,'TickLabelInterpreter','latex')

nexttile;
hold on

for i = 1:num_corners
    plot(t,u_H2(i,:),'linewidth',2)
end

title("Control Input (6MP Dose in mg)",'Interpreter','latex','FontSize',16)
xlabel('Time (days)','Interpreter','latex')
set(gca,'TickLabelInterpreter','latex')

P = value(X); % Lyapunov function parameter
save("robust_linear.mat", "pmax", "pmin", "K_H2", "P", "xbar_ave", "ubar_ave")
