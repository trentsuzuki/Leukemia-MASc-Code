%% Setup
clear;
clc;
close all;

if(~isdeployed)
  cd(fileparts(matlab.desktop.editor.getActiveFilename));
end

%% Simulate Adaptive Backstepping Control Law

x_d = 1.5;
k_ma = 2.5/3;

Base = 2.5;
k_tr = 0.15/3;
gamma_ = 0.8;
slope = 0.3/50;

Base_est = Base*1;
k_tr_est = k_tr*1.5;
slope_est = slope*1.5;

% Convergence/learning rates
c1 = k_ma / 2;
c2 = 1;
g1 = 0.00001;
g2 = 0.0001;
g3 = 0.0001;

% Projection function params
q = 6;
eps = 0.1;

dt = 0.01;
Tmax = 300;

t_values = 0:dt:Tmax;
max_idx = length(t_values);

convert_params = @(Base, k_tr, g, slope) [k_tr*Base^g; k_tr*slope*Base^g; k_tr];

theta_min = [0.001; 0.0001; 0.001];
theta_max = [1; 0.01; 1];

x = nan(max_idx,2);
z = nan(max_idx,2);
theta_hat = nan(max_idx,3);

theta_hat(1,:) = convert_params(Base_est, k_tr_est, gamma_, slope_est);
disp('theta_hat(1,:):')
disp(theta_hat(1,:))

theta = convert_params(Base, k_tr, gamma_, slope);


x(1,:) = [Base, k_tr/k_ma*Base];
% x(1,:) = [5, 0.09];
% x(1,:) = [5, 1];

alpha = nan(max_idx,1);
x_dot = nan(max_idx-1,2);
dalpha_dzeta1 = nan(max_idx-1,1);
dalpha_dtheta3_hat = nan(max_idx-1,1);

u = nan(max_idx,1);

u_eq = 1/slope*(1-Base^-gamma_*x_d^gamma_)
x2_eq = k_tr/k_ma*x_d

%% Simulation
    
for t = 1:length(t_values)-1
    alpha(t) = 1/theta_hat(t,3) * ((k_ma-c1)*x(t,1) + c1*x_d);

    z(t,1) = x(t,1) - x_d;
    z(t,2) = x(t,2) - alpha(t);

    dalpha_dzeta1(t) = 1/theta_hat(t,3) * (-c1 + k_ma);
    dalpha_dtheta3_hat(t) = -1/theta_hat(t,3)^2 * (-c1*z(t,1) + k_ma*x(t,1));
    
    theta1_hat_dot = g1*z(t,2)*x(t,2)*x(t,1)^-gamma_;
    theta3_hat_dot = g2 * ( z(t,1)*z(t,2) + z(t,1)/theta_hat(t,3) * (-c1*z(t,1) + k_ma*x(t,1)) - (1 + dalpha_dzeta1(t))*z(t,2)*x(t,2) );
    
    u(t) = 1/theta_hat(t,2)*( theta_hat(t,1) - theta_hat(t,3)*x(t,1)^gamma_ * (1 + dalpha_dzeta1(t)) + x(t,1)^gamma_/x(t,2) * (c2*z(t,2) + k_ma*dalpha_dzeta1(t)*x(t,1) - dalpha_dtheta3_hat(t)*theta3_hat_dot + theta_hat(t,3)*z(t,1)) );
    u(t) = min(max(u(t), 0), 300);% + sin(0.001*t)+1;
    
    theta2_hat_dot = -g3*u(t)*z(t,2)*x(t,2)*x(t,1)^-gamma_;
    
    % Project update laws to ensure parameters remain within bounds
    theta_hat_dot_proj = proj(theta_hat(t,:)', [theta1_hat_dot; theta2_hat_dot; theta3_hat_dot], theta_min, theta_max, eps, q);
    theta1_hat_dot = theta_hat_dot_proj(1);
    theta2_hat_dot = theta_hat_dot_proj(2);
    theta3_hat_dot = theta_hat_dot_proj(3);

    theta_hat(t+1,1) = theta_hat(t,1) + theta1_hat_dot*dt;
    theta_hat(t+1,2) = theta_hat(t,2) + theta2_hat_dot*dt;
    theta_hat(t+1,3) = theta_hat(t,3) + theta3_hat_dot*dt;

    x_dot(t,1) = theta(3)*x(t,2) - k_ma*x(t,1);
    x_dot(t,2) = theta(1)*x(t,1)^-gamma_*x(t,2) - theta(2)*x(t,1)^-gamma_*x(t,2)*u(t) - theta(3)*x(t,2);
    x(t+1,:) = x(t,:) + x_dot(t,:).*dt;
end

alpha(t+1) = 1/theta_hat(t+1,3) * ((k_ma-c1)*x(t+1,1) + c1*x_d);
z(t+1,1) = x(t+1,1) - x_d;
z(t+1,2) = x(t+1,2) - alpha(t+1);

%% Plot

f = figure('Position',[100 0 600 600]);
tiledlayout(8,1,"TileSpacing","tight","Padding","compact");

nexttile;
plot(t_values,x(:,1))
ylabel('$x_{ma}$','Interpreter','latex','Rotation',0,'FontSize',16)
nexttile;
plot(t_values,x(:,2))
ylabel('$x_{tr}$','Interpreter','latex','Rotation',0,'FontSize',16)
nexttile;
plot(t_values,z(:,1))
ylabel('$z_1$','Interpreter','latex','Rotation',0,'FontSize',16)
nexttile;
plot(t_values,z(:,2))
ylabel('$z_2$','Interpreter','latex','Rotation',0,'FontSize',16)
nexttile;
plot(t_values,theta_hat(:,1))
ylabel('$\hat\theta_1$','Interpreter','latex','Rotation',0,'FontSize',16)
nexttile;
plot(t_values,theta_hat(:,2))
ylabel('$\hat\theta_2$','Interpreter','latex','Rotation',0,'FontSize',16)
nexttile;
plot(t_values,theta_hat(:,3))
ylabel('$\hat\theta_3$','Interpreter','latex','Rotation',0,'FontSize',16)
nexttile;
plot(t_values,u)
ylabel('$u$','Interpreter','latex','Rotation',0,'FontSize',16)
xlabel('Time (days)')

current_time = string(datetime,'yyyyMMdd_HHmm');
% exportgraphics(f, 'simResult_' + current_time + '.png')
% close(f)

figure;
plot(x(:,1),x(:,2))
title('Phase portrait')
xlabel('x_1')
ylabel('x_2')

figure;
theta_tilde = theta_hat - theta';
V = 1/2*(z(:,1).^2 + z(:,2).^2 + 1/g1*theta_tilde(:,1).^2 + 1/g2*theta_tilde(:,2).^2 + 1/g3*theta_tilde(:,2).^2);
plot(t_values,V,'x')
title('Lyapunov Function')
xlabel('t')
ylabel('V')

%% Helper functions

function P = convex_param_set(p, pmin, pmax, eps, q)
    rho = (pmax + pmin) / 2;
    sigma = (pmax - pmin) / 2;
    terms = (p - rho) ./ sigma;
    P = 2/eps * (sum(terms.^q) - 1 + eps);
end

function tau_proj = proj(theta_hat, tau, theta_min, theta_max, eps, q)
    % Continuous approximation of projection function
    rho = (theta_max + theta_min) / 2;
    sigma = (theta_max - theta_min) / 2;
    P = convex_param_set(theta_hat,theta_min,theta_max,eps,q);
    dPdp = 2/eps* q * ((theta_hat - rho) ./ sigma).^(q-1)';
    
    tau_proj = tau;
    if P > 0 && dPdp*tau > 0
        tau_proj = tau_proj - (P*dPdp*tau)/(dPdp*dPdp') * dPdp';
    end
end

function tau_proj = proj_alt(theta_hat, tau, theta_min, theta_max, ~, ~)
    % Component-wise projection of updates to enforce box constraints
    tau_proj = tau;
    for i = 1:length(theta_hat)
        if theta_hat(i) <= theta_min(i) && tau_proj(i) < 0
            tau_proj(i) = 0;
        elseif theta_hat(i) >= theta_max(i) && tau_proj(i) > 0
            tau_proj(i) = 0;
        end
    end
end

function theta_proj = project_theta(theta_hat, theta_min, theta_max)
    % Component-wise projection of parameters to enforce box constraints
    theta_proj = theta_hat;
    for i = 1:length(theta_hat)
        if theta_hat(i) <= theta_min(i)
            theta_proj(i) = theta_min(i);
        elseif theta_hat(i) >= theta_max(i)
            theta_proj(i) = theta_max(i);
        end
    end
end
