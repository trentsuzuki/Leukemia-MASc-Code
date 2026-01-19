%% Setup
clear;
clc;

if(~isdeployed)
  cd(fileparts(matlab.desktop.editor.getActiveFilename));
end

%% Simulate Adaptive Backstepping Control Law

x_d = 1.5;
k_ma = 2.5;

Base = 2.5;
k_tr = 0.15;
gamma_ = 0.8;
slope = 0.3;

Base_est = 3;
k_tr_est = 0.2;
slope_est = 0.2;

c1 = 1;
c2 = 0.1;
mu = 1;
rho = 1;

dt = 0.001;
Tmax = 3000;

t_values = 0:dt:Tmax;
max_idx = length(t_values);

x = nan(max_idx,2);
x(1,1) = 2;
x(1,2) = 20;

x_dot = nan(max_idx-1,2);

u_eq = 1/slope*(1-Base^-gamma_*x_d^gamma_)
u = u_eq*ones(max_idx,1);

%% 
    
for t = 1:length(t_values)-1
    x_dot(t,1) = k_tr*x(t,2) - k_ma*x(t,1);
    x_dot(t,2) = k_tr*(1-slope*u(t))*Base^gamma_*x(t,1)^-gamma_*x(t,2) - k_tr*x(t,2);
    x(t+1,:) = x(t,:) + x_dot(t,:).*dt;
end

%% Plot

f = figure('Position',[1000 100 600 1200]);
tiledlayout(2,1,"TileSpacing","tight","Padding","compact");

nexttile;
plot(t_values,x(:,1))
ylabel('$x_{ma}$','Interpreter','latex','Rotation',0,'FontSize',16)
nexttile;
plot(t_values,x(:,2))
ylabel('$x_{tr}$','Interpreter','latex','Rotation',0,'FontSize',16)
xlabel('Time (days)')

current_time = string(datetime,'yyyyMMdd_HHmm');
% exportgraphics(f, 'simResult_' + current_time + '.png')
% close(f)