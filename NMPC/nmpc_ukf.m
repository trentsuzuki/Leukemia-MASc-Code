%% Nonlinear model predictive control using Unscented Kalman Filter
% Simplified pharmacokinetics compared to Jost 2020

clear; clc; close all;

%% Model/PK parameters

p.k_tr  = 0.15; % 1/day
p.k_ma  = 2.5; % 1/day
p.Base  = 2.5; % G/L
p.gamma = 0.8;
p.SLOPE = 0.27; % cells^-1*L
p.ke    = 0.2; % 1/day elimination of 6-MP surrogate concentration
p.F_over_V = 0.2; % lumped oral bioavailability / volume term

%% Discrete-time plant for MPC (sample = 7 days)

Ts_wk = 7; % days
nx = 6;
ny = 1;
nu = 1;
r = 1.5;

stateFcn  = @(x,u) log( friberg_integrator(exp(x),u,p,Ts_wk) );
outputFcn = @(x,u) exp(x(end)); % measurement is ANC only

%% 3. Build NL-MPC object

mpcobj = nlmpc(nx,ny,nu);
mpcobj.Ts = Ts_wk;
mpcobj.Model.StateFcn = stateFcn;
mpcobj.Model.OutputFcn = outputFcn;
mpcobj.Model.IsContinuousTime = false;
mpcobj.PredictionHorizon = 8; % weeks
mpcobj.ControlHorizon = 8; % weeks

% Constraints
mpcobj.MV.Min = 0; % mg/m^2/d
mpcobj.MV.Max = 75;
mpcobj.States(end).Min = 0.5; % safe lower bound for ANC
mpcobj.States(end).Max = 3; % safe upper bound for ANC

% Cost function weights
mpcobj.Weights.OutputVariables = 1; % Weight for ANC
mpcobj.Weights.ManipulatedVariables = 0.01; % Weight for drug dose
mpcobj.Weights.ManipulatedVariablesRate = 0; % Weight for drug dose
mpcobj.Weights.ECR = 1e5; % Slack variable tuning weight (default: 1e5)

validateFcns(mpcobj, rand(nx,1), 50);

% UKF on log transformed state to make noise additive
ukf = unscentedKalmanFilter( ...
    stateFcn, ...
    @(x) exp(x(end)), ...
    'HasAdditiveMeasurementNoise', true);

ukf.Alpha = 1e-2; % default 1e-3 can be too tight
ukf.Beta = 2; % "for Gaussian distributions, beta = 2 is optimal"
ukf.Kappa = 0; % "usually set to zero"

measurement_noise_sd = 0.1;
ukf.MeasurementNoise = measurement_noise_sd.^2; % Variance of ANC measurement noise (G/L)
ukf.ProcessNoise = 1e-5*eye(nx);

% Initialization
x0 = [0.1; p.Base*p.k_ma/p.k_tr*ones(4,1); p.Base];
ukf.State = log(x0);
ukf.StateCovariance = 0.01^2 * eye(nx);

%% 5. Closed-loop simulation

Nk = 50; % weeks to simulate
xtrue = x0;
dose_prev = 0; % initial daily dose (mg/m^2/d)

hist.y = zeros(1,Nk);
hist.u = zeros(1,Nk);
hist.t = Ts_wk*(0:Nk-1);

rng(0);
startTime = datetime('now');
for k = 1:Nk
    hist.x(k,:) = xtrue;

    % Measurement (clinic visit)
    yk = xtrue(end) * (1 + measurement_noise_sd*randn);
    hist.y(k) = yk;

    % UKF update with 7-day batch
    correct(ukf, yk);
    hist.xhat(k,:) = ukf.State;

    % NMPC optimisation
    [uk,info] = nlmpcmove(mpcobj, ukf.State, dose_prev, r);

    % Store & apply weekly constant dose
    hist.u(k) = uk;

    % Simulate true dynamics over next week
    dt = 1/24;
    Nt = Ts_wk/dt;
    for i = 1:Nt
        dose = 0;
        % Daily dosing event at hour 0 of each day
        timeofday = mod((i-1)*dt,1);
        if abs(timeofday) == 0
            dose = uk;
        end
        xtrue = rk4(@(z) friberg_rhs(z,dose,p), xtrue, dt);
    end

    % UKF time-update for the whole week (Ts_wk) with same uk
    predict(ukf, uk);

    dose_prev = uk;
end
stopTime = datetime('now');

disp(stopTime - startTime);

%% Plots

figure;
subplot(2,1,1);
hold on;
stairs(hist.t, hist.u, 'LineWidth',1.4);
xlabel('Day');
ylabel('6-MP dose (mg/m^2/d)');
grid on;
subplot(2,1,2);
hold on;
plot(hist.t, hist.y, 'r*');
stairs(hist.t, hist.x(:,6), 'LineWidth',1.4);
yline(mpcobj.States(5).Min,'r--','Lower safety bound');
yline(mpcobj.States(5).Max,'r--','Upper safety bound');
ylim([0 3])
xlabel('Day');
ylabel('ANC (10^9 cells/L)');
legend({'Measured','True'})
grid on;
title('Measured output')

figure;
hF = gcf;
hF.Position = [300 300 600 700];
subplot(3,1,1);
hold on;
stairs(hist.t, hist.u, 'LineWidth',1.4);
xlabel('Day');
ylabel('6-MP dose (mg/m^2/d)');
grid on;
title('Drug Dose')
subplot(3,1,2);
hold on;
stairs(hist.t, hist.x, 'LineWidth',1.4);
yline(mpcobj.States(5).Min,'r--','Lower safety bound');
yline(mpcobj.States(5).Max,'r--','Upper safety bound');
xlabel('Day');
ylabel('ANC (10^9 cells/L)');
grid on;
title('True states')
subplot(3,1,3);
hold on;
stairs(hist.t, exp(hist.xhat), 'LineWidth',1.4);
yline(mpcobj.States(5).Min,'r--','Lower safety bound');
yline(mpcobj.States(5).Max,'r--','Upper safety bound');
xlabel('Day');
ylabel('ANC (10^9 cells/L)');
grid on;
title('Estimated states')

%% Helper functions

% Integrated Friberg dynamics over the time period Delta
function xnext = friberg_integrator(x,u,p,Delta)
    % Integrate PK+PD over Delta (days) with daily bolus dose u.
    dt = 1/24; % internal step of 1 hour
    Nt = round(Delta/dt);
    for n = 1:Nt
        dose = 0;
        % Dose event happens once every day at hour zero
        if mod((n-1)*dt,1) == 0
            dose = u;
        end
        x = rk4(@(z) friberg_rhs(z,dose,p), x, dt);
    end
    xnext = x;
end

% PK/PK dynamics
function dx = friberg_rhs(x,dose,p)
    x = max(x,1e-8); % Ensure positive states

    C=x(1); P=x(2); T1=x(3); T2=x(4); T3=x(5); N=x(6);

    Edrug = p.SLOPE*C;
    feedback = (p.Base/N)^p.gamma;
    
    % Single compartment PK:
    dC  = -p.ke*C + p.F_over_V*dose;

    % Friberg model with 3 transitional states:
    dP  = p.k_tr*(1-Edrug)*P*feedback - p.k_tr*P;
    dT1 = p.k_tr*P  - p.k_tr*T1;
    dT2 = p.k_tr*T1 - p.k_tr*T2;
    dT3 = p.k_tr*T2 - p.k_tr*T3;
    dN  = p.k_tr*T3 - p.k_ma*N;

    dx = [dC; dP; dT1; dT2; dT3; dN];
end

% 4th order Runge-Kutta
% https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods
function xnew = rk4(f,x,h)
    k1 = f(x);
    k2 = f(x + 0.5*h*k1);
    k3 = f(x + 0.5*h*k2);
    k4 = f(x + h*k3);
    xnew = x + h/6*(k1 + 2*k2 + 2*k3 + k4);
end
