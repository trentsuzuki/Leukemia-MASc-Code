% Simulation of Adaptive Backstepping for the simplified Jost model
% with Parameter Projection to fix the controller singularity flaw.

clear; clc;
close all;

set(0,'defaultTextInterpreter','latex');

%% Define Known Parameters and True (Unknown) Parameters

% Known
p.kma = 2.375;    % Death rate of mature cells (1/day) - Example value
p.gamma = 0.7;  % Feedback sensitivity - Example value

% True params (unknown by the controller)
p.k_tr = 0.15/2.8;
p.Base = 2.5;
p.slope = 0.3;

% Convert from Base, k_tr, slope to theta
convert_params = @(Base, k_tr, g, slope) [k_tr*Base^g; k_tr*slope*Base^g; k_tr];

theta_true = convert_params(p.Base, p.k_tr, p.gamma, p.slope);
p.theta1_true = theta_true(1);
p.theta2_true = theta_true(2);
p.theta3_true = theta_true(3);

% Controller Design Parameters
p.Nd = 1.5; % Desired setpoint for N

% Gains c1, c2. Note: c1 must satisfy 0 < c1 < kma
p.c1 = 0.8;
p.c2 = 1;

% Adaptation gains (Gamma matrix diagonal)
p.g1 = 1.0e-6;
p.g2 = 1.0e-7;
p.g3 = 1.0e-8;

% Projection Bounds
% These define the known convex set containing the true parameters.
Base_l = 1.5;
k_tr_l = 0.01;
slope_l = 0.2;
% Lower bounds MUST be strictly positive to avoid singularity.
Base_u = 6;
k_tr_u = 0.1;
slope_u = 0.6;

p.theta_min = convert_params(Base_l, k_tr_l, p.gamma, slope_l);
p.theta_max = convert_params(Base_u, k_tr_u, p.gamma, slope_u);

if any(p.theta_min <= 0)
    error('Lower projection bounds must be strictly positive.');
end

%% Simulation Setup

% Initial conditions [N, P, hat_theta1, hat_theta2, hat_theta3]
N_0 = p.Base; % Initial N (Baseline equilibrium condition)
P_0 = p.kma/p.k_tr*p.Base; % Initial P
% Initial estimates (must be within bounds)
theta_hat_0 = (p.theta_max + p.theta_min)/2; % Average of bounds

if any(theta_hat_0 < p.theta_min) || any(theta_hat_0 > p.theta_max)
    error('Initial estimates are outside the projection bounds.');
end

Y0 = [N_0; P_0; theta_hat_0];
tspan = [0 100]; % Simulation time

%% Run Simulation

disp('Starting simulation...');

options = odeset('RelTol',1e-6,'AbsTol',1e-6);
[t, Y] = ode45(@(t, y) closedLoopDynamics(t, y, p), tspan, Y0, options);

disp('Simulation finished.');

%% Extract Results and Visualize

N = Y(:,1);
P = Y(:,2);
theta_hat1 = Y(:,3);
theta_hat2 = Y(:,4);
theta_hat3 = Y(:,5);

% Recalculate control input u for plotting
u_plot = zeros(size(t));
for i = 1:length(t)
    [~, u_plot(i)] = closedLoopDynamics(t(i), Y(i,:)', p);
end

plotResults(t, N, P, u_plot, theta_hat1, theta_hat2, theta_hat3, p);

%% Helper Functions

function [dYdt, u] = closedLoopDynamics(~, Y, p)
    % Implements the closed-loop dynamics (System + Controller + Adaptation)

    % Extract states and parameters
    N = Y(1);
    P = Y(2);
    theta_hat = Y(3:5);
    theta_hat1 = theta_hat(1);
    theta_hat2 = theta_hat(2);
    theta_hat3 = theta_hat(3);

    kma = p.kma;
    gamma = p.gamma;
    Nd = p.Nd;
    c1 = p.c1;
    c2 = p.c2;

    % Safety check for state positivity (required for N^-gamma)
    % Prevents numerical issues if states approach zero.
    if N <= 1e-8
        N = 1e-8;
    end
     if P <= 1e-8
        P = 1e-8;
    end

    % --- Control Design Implementation ---

    % Error z1
    z1 = N - Nd;

    % Virtual control alpha and its partial derivatives
    % alpha = (1/h_th3) * ((kma-c1)*N + c1*Nd)
    alpha = (1/theta_hat3) * ((kma - c1)*N + c1*Nd);

    % d_alpha/d_N
    d_alpha_d_N = (kma - c1) / theta_hat3;

    % d_alpha/d_h_th3
    d_alpha_d_th3 = -alpha / theta_hat3;

    % Error z2 (z2 = exp(xi2) - alpha = P - alpha)
    z2 = P - alpha;

    % Calculate Tuning Function tau_3 (Independent of u)
    term1 = z1 * z2;
    term2 = (z1/theta_hat3)*(-c1*z1 + kma*N);
    % The term below combines the contribution from phi(3)=-1 and the backstepping derivative.
    term3 = -(1 + d_alpha_d_N)*z2*P;

    tau3 = p.g3 * (term1 + term2 + term3);

    % Calculate d_theta_hat3 using Projection
    % This MUST be done before calculating u, as u depends on it.
    d_theta_hat3 = projectionOperator(tau3, theta_hat3, p.theta_min(3), p.theta_max(3));

    % Calculate Control Law u
    N_g = N^gamma;

    U_term1 = theta_hat1 - theta_hat3*N_g*(1 + d_alpha_d_N);
    % Note: We use the projected d_theta_hat3 here
    U_term2_inner = c2*z2 + kma*d_alpha_d_N*N - d_alpha_d_th3*d_theta_hat3 + theta_hat3*z1;
    U_term2 = (N_g / P) * U_term2_inner;

    % Division by theta_hat2 is now safe due to projection
    u = (1/theta_hat2) * (U_term1 + U_term2);

    % Physical constraint: Drug concentration cannot be negative.
    % u = max(0, u);

    % Calculate Tuning Functions tau_1 and tau_2 (Dependent on u)
    N_neg_g = 1/N_g;
    common_term = z2 * P * N_neg_g;

    tau1 = p.g1 * common_term;
    tau2 = -p.g2 * common_term * u;

    % Calculate d_theta_hat1 and d_theta_hat2 using Projection
    d_theta_hat1 = projectionOperator(tau1, theta_hat1, p.theta_min(1), p.theta_max(1));
    d_theta_hat2 = projectionOperator(tau2, theta_hat2, p.theta_min(2), p.theta_max(2));

    d_theta_hat = [d_theta_hat1; d_theta_hat2; d_theta_hat3];

    % --- System Dynamics (for simulation using true parameters) ---
    theta1 = p.theta1_true;
    theta2 = p.theta2_true;
    theta3 = p.theta3_true;

    % dN = theta3*P - kma*N
    dNdt = theta3*P - kma*N;
    % dP = theta1*N^(-gamma)*P - theta2*N^(-gamma)*P*u - theta3*P
    dPdt = theta1*N_neg_g*P - theta2*N_neg_g*P*u - theta3*P;

    dYdt = [dNdt; dPdt; d_theta_hat];
end

function d_theta_i = projectionOperator(tau_i, theta_i, min_val, max_val)
    % Implements the projection algorithm for scalar parameters
    % Ensures the parameter estimate remains within [min_val, max_val]
    if (theta_i <= min_val && tau_i < 0)
        % At lower bound and trying to decrease
        d_theta_i = 0;
    elseif (theta_i >= max_val && tau_i > 0)
        % At upper bound and trying to increase
        d_theta_i = 0;
    else
        % Within bounds or moving inwards
        d_theta_i = tau_i;
    end
end

function plotResults(t, N, P, u_plot, theta_hat1, theta_hat2, theta_hat3, p)
    % Helper function for plotting the results
    figure('Name', 'Adaptive Backstepping with Projection', 'Position', [100, 100, 1000, 800], 'Theme', 'Light');

    % Plot 1: System output
    subplot(3,1,1);
    hold on;
    plot(t, N, 'r', 'LineWidth', 2);
    % plot(t, P, 'b--', 'LineWidth', 2);
    yline(p.Nd, 'k--', 'LineWidth', 1.5);
    title('System States');
    xlabel('Time (days)');
    ylabel('Concentration');
    legend('N (Mature Neutrophils)', ... 'x_2 (Proliferating Cells)', ...
        'y_r (Reference)', ...
        'Location', 'best');
    grid on;

    % Plot 2: Control Input
    subplot(3,1,2);
    plot(t, u_plot, 'm', 'LineWidth', 2);
    plot(t, u_plot*p.slope, 'm', 'LineWidth', 2);
    title('Control Input (u)');
    xlabel('Time (days)');
    % ylabel('u (6-TGN Concentration)');
    ylabel('u (drug effect)');
    ylim([0 ceil(max(u_plot*p.slope))])
    grid on;

    % Plot 3: Parameter Estimates
    subplot(3,1,3);
    hold on;
    plot(t, theta_hat1, 'r', 'LineWidth', 2);
    plot(t, theta_hat2, 'g', 'LineWidth', 2);
    plot(t, theta_hat3, 'b', 'LineWidth', 2);
    % Plot true values with dotted lines
    yline(p.theta1_true, 'r:', 'LineWidth', 1.5);
    yline(p.theta2_true, 'g:', 'LineWidth', 1.5);
    yline(p.theta3_true, 'b:', 'LineWidth', 1.5);
    title('Parameter Estimates ($\hat\theta$)', 'Interpreter', 'latex');
    xlabel('Time (days)');
    ylabel('Parameter Estimate');
    legend('$\hat\theta_1$', '$\hat\theta_2$', '$\hat\theta_3$', '$\theta$', 'Location', 'east', 'Interpreter','Latex');
    grid on;
end
