%% Robust control of the model from Jost 2020
% State feedback using H infinity synthesis

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
pmin.slope = 0.25;
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

%% Rescaling
% Compute scaling factors based on the equilibrium state xbar_ave.
% Use a small threshold (e.g. 1e-3) in case any entry is nearly zero.
% scale_factors = max(abs(xbar_ave), 1e-2);  
scale_factors = [46.4172   87.7086    3.3486   1/k_tr_sample 1/k_tr_sample 1/k_tr_sample 1/k_tr_sample   0.9193];
% scale_factors = ones(n,1);
T = diag(1./scale_factors);  % Scaling matrix, so that scaled states are ~1 in magnitude
% T = diag(1);  % Scaling matrix, so that scaled states are ~1 in magnitude

disp('Scaling factors for each state:');
disp(scale_factors);

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
A_scaled = cell(num_corners,1);

% Generate all corners of the convex set
for i = 0:num_corners-1
    binary_index = dec2bin(i, num_dimensions) - '0'; % Convert to binary array
    corners(i+1, :) = delta_min .* (1 - binary_index) + delta_max .* binary_index;
    A_scaled{i+1} = T * A(corners(i+1,:)) * inv(T);
    eig(A_scaled{i+1});
end

Bu = [F; zeros(7,1)];
Cy = [zeros(1,7) 1];

sys_nom = compreal(ss(A_nom,Bu,Cy,[]),"o"); % observable canonical form

%% Robust H infinity synthesis

m = 1;
nz = 2;
nw = 1;

Bw = x0_jost - xbar_ave;
Cz = [zeros(1,n-1) 1e-3; zeros(m,8)];
Dzw = zeros(nz,nw);
Dzu = [zeros(1,m); 1e-3];

X = sdpvar(n,n);
Z = sdpvar(m,n,'full');
gamma_sq = sdpvar;
small = 5e-1;
eps = 1e-5;

Bu_scaled = T * Bu;
Cz_scaled = Cz * inv(T);

M = cell(num_corners,1);

% Generate all corners of the convex set
for i = 1:num_corners
    M{i} = [
        (A_scaled{i}*X+Bu_scaled*Z)' + A_scaled{i}*X+Bu_scaled*Z, Bw, (Cz_scaled*X+Dzu*Z)';
        Bw', -gamma_sq*eye(nw), Dzw';
        (Cz_scaled*X+Dzu*Z), Dzw, -eye(nz);
    ];
end

Constraints = [
    X >= small*eye(n);
    M{1} <= -eps*eye(n+nw+nz);
    M{2} <= -eps*eye(n+nw+nz);
    M{3} <= -eps*eye(n+nw+nz);
    M{4} <= -eps*eye(n+nw+nz);
    M{5} <= -eps*eye(n+nw+nz);
    M{6} <= -eps*eye(n+nw+nz);
    M{7} <= -eps*eye(n+nw+nz);
    M{8} <= -eps*eye(n+nw+nz);
    M{9} <= -eps*eye(n+nw+nz);
    M{10} <= -eps*eye(n+nw+nz);
    M{11} <= -eps*eye(n+nw+nz);
    M{12} <= -eps*eye(n+nw+nz);
    M{13} <= -eps*eye(n+nw+nz);
    M{14} <= -eps*eye(n+nw+nz);
    M{15} <= -eps*eye(n+nw+nz);
    M{16} <= -eps*eye(n+nw+nz);
    M{17} <= -eps*eye(n+nw+nz);
    M{18} <= -eps*eye(n+nw+nz);
    M{19} <= -eps*eye(n+nw+nz);
    M{20} <= -eps*eye(n+nw+nz);
    M{21} <= -eps*eye(n+nw+nz);
    M{22} <= -eps*eye(n+nw+nz);
    M{23} <= -eps*eye(n+nw+nz);
    M{24} <= -eps*eye(n+nw+nz);
    M{25} <= -eps*eye(n+nw+nz);
    M{26} <= -eps*eye(n+nw+nz);
    M{27} <= -eps*eye(n+nw+nz);
    M{28} <= -eps*eye(n+nw+nz);
    M{29} <= -eps*eye(n+nw+nz);
    M{30} <= -eps*eye(n+nw+nz);
    M{31} <= -eps*eye(n+nw+nz);
    M{32} <= -eps*eye(n+nw+nz);
];

Cost = gamma_sq;
options = sdpsettings('solver','sdpt3','verbose',1,'sdpt3.maxit',1000,'sdpt3.gam',0.6);
sol = optimize(Constraints,Cost,options);

sol.info
sqrt(value(Cost))
K_scaled = value(Z)*inv(value(X))
K_Hinf = K_scaled * T;

%% Simulate H infinity controller on linear system

t = 0:0.1:400;
w = zeros(length(t),nw);

x_Hinf = cell(num_corners,1);
u_Hinf = nan(num_corners,length(t));
output_names = {'x_1','x_2','x_3','x_4','x_5','x_6','x_7','x_8'};

for i = 1:num_corners
    sys_i = ss(A(corners(i,:))+Bu*K_Hinf,zeros(n,nw),eye(n),zeros(n,nw),'TimeUnit','Days','OutputName',output_names);
    x_Hinf(i) = {lsim(sys_i,w,t,x0_jost-xbar_ave)};
    u_Hinf(i,:) = K_Hinf*x_Hinf{i}' + ubar_ave;
end

f_Hinf = figure('Position',[0 0 800 500]);
tiledlayout(2,1,"TileSpacing","tight","Padding","compact");

nexttile;
hold on

for i = 1:num_corners
    x = x_Hinf{i};
    plot(t,x(:,8) + xbar_8,'linewidth',2)
end

yline(xbar_8,'k--')
title("Disturbance Impulse Response (Neutrophil Count in G/L)",'Interpreter','latex','FontSize',16)
xlabel('Time (days)','Interpreter','latex')
set(gca,'TickLabelInterpreter','latex')

nexttile;
hold on

for i = 1:num_corners
    plot(t,u_Hinf(i,:),'linewidth',2)
end

title("Control Input (6MP Dose in mg)",'Interpreter','latex','FontSize',16)
xlabel('Time (days)','Interpreter','latex')
set(gca,'TickLabelInterpreter','latex')

% P = value(X); % Lyapunov function parameter
% save("robust_Hinf.mat", "pmax", "pmin", "K_Hinf", "P", "xbar_ave", "ubar_ave")

for i = 1:num_corners
    disp(max(abs(x_Hinf{i})))
end