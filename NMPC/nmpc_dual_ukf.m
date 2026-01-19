%% Nonlinear model predictive control using Dual Unscented Kalman Filter
% Simplified pharmacokinetics compared to Jost 2020

clear; clc; close all;

%% Model/PK parameters

p.ke = 0.2; % 1/day elimination of 6-MP surrogate concentration
p.F_over_V = 0.2; % lumped oral bioavailability / volume term
p.Base = 2.5; % G/L
p.k_tr = 0.15; % 1/day
p.SLOPE = 0.27; % cells^-1*L
p.gamma = 0.8;
p.k_ma = 2.5; % 1/day

% Parameter estimates (only Base, k_tr, SLOPE, and gamma are adapted)
phat.ke = 0.2;
phat.F_over_V = 0.2;

paramEstFactor = 1.2;
phat.Base = 2.5*paramEstFactor; % G/L
phat.k_tr = 0.15*paramEstFactor; % 1/day
phat.SLOPE = 0.27*paramEstFactor; % cells^-1*L
phat.gamma = 0.8*paramEstFactor;

phat.k_ma = 2.5; % 1/day

%% Discrete-time plant for MPC (sample = 7 days)

Ts_wk = 7; % days
nx = 6;
ny = 1;
nu = 1;
np = 4;
r = 1.25; % target ANC

stateFcn = @(x,u) zeros(nx,1); % Placeholder to be updated later
mpcOutputFcn = @(x,u) exp(x(nx)); % Measurement is ANC only

%% Build NL-MPC object

mpcobj = nlmpc(nx,ny,nu);
mpcobj.Ts = Ts_wk;
mpcobj.Model.StateFcn = stateFcn; % Updated during sim to use the latest data
mpcobj.Model.OutputFcn = mpcOutputFcn;
mpcobj.Model.IsContinuousTime = false;
mpcobj.PredictionHorizon = 8; % weeks
mpcobj.ControlHorizon = 8; % weeks

% Constraints
mpcobj.MV.Min = 0; % mg/m^2/d
mpcobj.MV.Max = 75;
% mpcobj.States(nx).Min = 0.5; % safe lower bound for ANC
% mpcobj.States(nx).Max = 2; % safe upper bound for ANC

% Cost function weights
mpcobj.Weights.OutputVariables = 1; % Weight for ANC
mpcobj.Weights.ManipulatedVariables = 0.01; % Weight for drug dose
mpcobj.Weights.ManipulatedVariablesRate = 0; % Weight for drug dose
mpcobj.Weights.ECR = 1e5; % Slack variable tuning weight (default: 1e5)

validateFcns(mpcobj, rand(nx,1), 50);

%% Build Unscented Kalman Filter for state estimation

% UKF on log transformed plant state to make noise additive
% Assuming multiplicative measurement noise on parameters as well
stateUkf = unscentedKalmanFilter( ...
    stateFcn, ...
    @(x) exp(x(nx)), ...
    'HasAdditiveMeasurementNoise', true);

stateUkf.Alpha = 1e-2; % default 1e-3 can be too tight
stateUkf.Beta = 2; % "for Gaussian distributions, beta = 2 is optimal"
stateUkf.Kappa = 0; % "usually set to zero"

measurement_noise_sd = 0.1;
stateUkf.MeasurementNoise = measurement_noise_sd.^2; % Variance of ANC measurement noise (G/L)
stateUkf.ProcessNoise = 1e-5*eye(nx);

% Initialization
xhat0 = [0.1; phat.Base*phat.k_ma/phat.k_tr*ones(4,1); phat.Base];
stateUkf.State = log(xhat0);
stateUkf.StateCovariance = 1e-2^2*eye(nx);

%% Build Unscented Kalman Filter for parameter estimate

paramUkfOutputFcn = @(x) exp(x(nx));

% UKF on log transformed plant state to make noise additive
% Assuming multiplicative measurement noise on parameters as well
paramUkf = unscentedKalmanFilter( ...
    @(x,u) x, ...
    @(x) paramUkfOutputFcn(stateFcn(x,uk)), ...
    'HasAdditiveMeasurementNoise', true);

paramUkf.Alpha = 1e-2; % default 1e-3 can be too tight
paramUkf.Beta = 2; % "for Gaussian distributions, beta = 2 is optimal"
paramUkf.Kappa = 0; % "usually set to zero"

param_meas_noise_sd = 0.1;
paramUkf.MeasurementNoise = param_meas_noise_sd.^2; % Variance of ANC measurement noise (G/L)
paramUkf.ProcessNoise = 1e-5*eye(np);

% Initialization
phat0 = [phat.Base; phat.k_tr; phat.SLOPE; phat.gamma]; % phat
paramUkf.State = log(phat0);
paramUkf.StateCovariance = 1e-1^2*eye(np);

%% 5. Closed-loop simulation

Nk = 30; % Weeks to simulate
dose_prev = 0; % Initial daily dose (mg/m^2/d)
x0 = [0.1; p.Base*p.k_ma/p.k_tr*ones(4,1); p.Base]; % Initial state
xtrue = x0;

hist.x = nan(Nk,nx); % stored in linear scale
hist.xhat = nan(Nk,nx); % stored with log transform applied
hist.phat = nan(Nk,np); % stored with log transform applied
hist.y = nan(Nk,1); % stored in linear scale
hist.u = nan(Nk,1); % stored in linear scale
hist.t = Ts_wk*(0:Nk-1);

hist.u(1) = dose_prev;
hist.xhat(1,:) = log(xhat0);
hist.phat(1,:) = log(phat0);
hist.x(1,:) = x0;
hist.y(1,:) = hist.x(1,nx) * (1 + measurement_noise_sd*randn);

rngSeed = 1;
rng(rngSeed);
startTime = datetime('now');
for k = 2:Nk
    disp("Iteration: " + (k-1) + "/" + (Nk-1));
    
    hist.x(k,:) = xtrue;

    % Measurement (clinic visit)
    hist.y(k) = hist.x(k,nx) * (1 + measurement_noise_sd*randn);

    % State UKF update with 7-day batch:

    % Update state function to use the latest parameter estimate
    stateUkfUpdateFcn = @(x,u) log( friberg_integrator(exp(x),u,makeParamStruct(p,exp(hist.phat(k-1,:))),Ts_wk) );
    stateUkf.StateTransitionFcn = stateUkfUpdateFcn;

    correct(stateUkf, hist.y(k));
    hist.xhat(k,:) = stateUkf.State;
    
    % Parameter UKF update with 7-day batch:
    paramUkfMeasFcn = @(phat) log( friberg_integrator(exp(hist.xhat(k-1,:)'),hist.u(k-1),makeParamStruct(p,exp(phat)),Ts_wk) );
    paramUkf.MeasurementFcn = @(phat) paramUkfOutputFcn(paramUkfMeasFcn(phat));
    
    correct(paramUkf, hist.y(k));
    hist.phat(k,:) = paramUkf.State;

    % NMPC optimisation
    mpcStateFcn = @(x,u) log( friberg_integrator(exp(x),u,makeParamStruct(p,exp(hist.phat(k,:))),Ts_wk) );
    mpcobj.Model.StateFcn = mpcStateFcn; % Update state function to use latest estimates
    [uk,info] = nlmpcmove(mpcobj, hist.xhat(k,:), dose_prev, r);

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
        xtrue = rk4(@(x) friberg_rhs(x,dose,p), xtrue, dt);
    end

    % UKF time-update for the whole week (Ts_wk) with same uk
    predict(stateUkf, uk);
    predict(paramUkf, uk);

    dose_prev = uk;
end
stopTime = datetime('now');

disp("Simulation finished. Duration:");
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
stairs(hist.t, hist.x(:,nx), 'LineWidth',1.4);
yline(r,'b--','Target');
yline(mpcobj.States(nx).Min,'r--','Lower safety bound');
yline(mpcobj.States(nx).Max,'r--','Upper safety bound');
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
xlabel('Day');
grid on;
title('True states')
subplot(3,1,3);
hold on;
stairs(hist.t, exp(hist.xhat), 'LineWidth',1.4);
xlabel('Day');
grid on;
title('Estimated states')

figure;
hF = gcf;
hF.Position = [300 300 600 700];
hold on;
stairs(hist.t, exp(hist.phat), 'LineWidth',1.4);
yline(p.Base,'--r','Base')
yline(p.k_tr,'--r','k_{tr}')
yline(p.SLOPE,'--r','slope')
yline(p.gamma,'--r','gamma')
xlabel('Day');
grid on;
title('Parameter Estimates')

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

function pstruct = makeParamStruct(p,pvec)
    pstruct.ke = p.ke;
    pstruct.F_over_V = p.F_over_V;
    pstruct.Base = pvec(1);
    pstruct.k_tr = pvec(2);
    pstruct.SLOPE = pvec(3);
    pstruct.gamma = pvec(4);
    pstruct.k_ma = p.k_ma;
end
