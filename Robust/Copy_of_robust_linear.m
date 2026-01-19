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

%% ==========  Uncertainty description with algebraic coupling  ==========
%
%  Independent patient parameters
%      theta1 = Base      (G/L)
%      theta2 = k_tr      (1 / day)
%      theta3 = slope     (–)
%      theta4 = γ         (–)
%
%  All other factors in A(theta) are computed from these four.
%

% ----------  bounds and nominal values  ----------
thmin = [pmin.Base ; pmin.k_tr ; pmin.slope ; pmin.gamma];
thmax = [pmax.Base ; pmax.k_tr ; pmax.slope ; pmax.gamma];
thnom = 0.5*(thmin+thmax);           % mid-point – used only for plots etc.

% ----------  helper: coupled Delta-terms ----------
% Implements the algebraic relations once, so we can re-use it everywhere
theta2delta = @(th) [ ...
       th(2);                                                      % Delta1 = k_tr
       th(2)*(th(1)/xbar_8)^th(4);                                 % Delta2
       th(2)*th(4)*(th(1)/xbar_8)^th(4);                           % Delta3
       th(2)*th(4)*(th(1)/xbar_8)^th(4)*th(3);                     % Delta4
       th(2)*(th(1)/xbar_8)^th(4)*th(3) ];                         % Delta5

% ----------  A(theta) still affine – but only in the *independent* thetai ----------
fi  = @(th) theta2delta(th);       % shorthand
Atheta  = @(th) A0 + [A1,A2,A3,A4,A5]*fi(th);

% ----------  LFT realisation ------------------------------------------------
% Make a (static) LFT  M + N*Delta*L  where |Delta_i| ≤ 1.
% 1. Normalise every thetai to the box [-1,1]
Stheta  = diag( 2./(thmax-thmin) );    % scaling
thetabar = -ones(4,1) + Stheta*( thnom - thmin );   % normalised nominal (≈0)

% 2. Linear fractional transform
Delta   = ultidyn('Delta', [4 4]);
Delta.SamplingGrid = struct;     % diag{Delta_i}
Delta   = blkdiag(Delta(:,:,1),Delta(:,:,2),Delta(:,:,3),Delta(:,:,4));          % static, ‖Delta‖∞≤1

theta   = thnom + 0.5*(thmax-thmin).*Delta;      % box-to-ball mapping

% 3. *All* non-affine factors are now algebraic functions of theta,
%    so A(theta) can be written as an lft:
A_LFT = lft(Delta, Atheta(theta));          % returns a ‘uss’ object (nominal+Delta)

% ----------  Generate the 2^4 = 16 real vertices for the LMIs  ----------
num_par     = 4;
num_corners = 2^num_par;
A_corners   = cell(num_corners,1);
corners     = zeros(num_corners,num_par);

for i = 0:num_corners-1
    idx            = dec2bin(i,num_par)-'0';        % 0/1 for min/max
    th_i           = thmin.*(~idx) + thmax.*idx;    % actual theta-corner
    corners(i+1,:) = th_i.';
    A_corners{i+1} = Atheta(th_i);                      % coupled A-matrix
end

% Optional: sanity-check eigenvalues
% cellfun(@(Ai) max(real(eig(Ai))),A_corners)
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
