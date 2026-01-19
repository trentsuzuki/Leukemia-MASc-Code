%% Nonlinear model predictive control using Unscented Kalman Filter (No Parameter Adaptation)
% Simplified pharmacokinetics compared to Jost 2020

clear; clc; close all;

% Path to CSV file with patient parameters (Data Sheet 1 from Jost 2020 supplementary materials)
csvFilePath = "C:\Users\trent\UofT\Research\Leukemia\Code\Jost 2020 Data Analysis\Supplementary Materials\Data Sheet 1.CSV";

%% Read patient parameters from CSV
patientData = readtable(csvFilePath);
numPatients = height(patientData);

fprintf('Found %d patients in CSV file\n', numPatients);

%% Parallel pool setup
% Check if parallel pool is already running
poolobj = gcp('nocreate');
if isempty(poolobj)
    fprintf('Starting parallel pool...\n');
    parpool('local'); % Start local parallel pool
    poolobj = gcp;
end
fprintf('Using parallel pool with %d workers\n', poolobj.NumWorkers);

%% Fixed model/PK parameters (same for all patients)
p_fixed.ke = 0.2; % 1/day elimination of 6-MP surrogate concentration
p_fixed.F_over_V = 0.2; % lumped oral bioavailability / volume term
p_fixed.k_ma = 2.5; % 1/day

ancMeasNoiseSD = 0.1;

%% Simulation parameters
Ts_wk = 7; % days
nx = 6;
ny = 1;
nu = 1;
np = 4;
Nk = 30; % weeks to simulate

ancSafeLowerLim = 0.5; % If ANC is below this level, give zero dose

% MPC parameters
mpcAncTarget = 1.5; % target ANC
mpcMaxAllowableDose = 100; % Maximum 6-MP dose (mg/m^2/d)
mpcAncWeight = 1; % Weight for ANC
mpcDoseWeight = 0.01; % Weight for 6-MP dose
mpcHorizon = 8; % Number of weeks to predict and control
mpcECR = 1e5; % Slack variable tuning weight

% UKF parameters (Important: these are in the log domain)
ukfStateProcessNoiseCovs = [1e-6 1e-4 1e-4 1e-4 1e-4 1e-6];
ukfInitialStateSDs = [1e-4 0.125304565*ones(1,4) 0.103687761];

%% Validate functions once before parallel execution
fprintf('Validating MPC functions...\n');
p_temp = p_fixed;
p_temp.Base = 2.5;
p_temp.k_tr = 0.15;
p_temp.slope = 0.27;
p_temp.gamma = 0.8;

stateFcn_temp = @(z,u) log( friberg_integrator(exp(z),u,p_temp,Ts_wk) );
outputFcn_temp = @(x,u) exp(x(nx));

mpcobj_temp = nlmpc(nx,ny,nu);
mpcobj_temp.Ts = Ts_wk;
mpcobj_temp.Model.StateFcn = stateFcn_temp;
mpcobj_temp.Model.OutputFcn = outputFcn_temp;
mpcobj_temp.Model.IsContinuousTime = false;
mpcobj_temp.PredictionHorizon = mpcHorizon;
mpcobj_temp.ControlHorizon = mpcHorizon;
mpcobj_temp.MV.Min = 0;
mpcobj_temp.MV.Max = mpcMaxAllowableDose;
mpcobj_temp.Weights.OutputVariables = mpcAncWeight;
mpcobj_temp.Weights.ManipulatedVariables = mpcDoseWeight;
mpcobj_temp.Weights.ManipulatedVariablesRate = 0;
mpcobj_temp.Weights.ECR = 1e5;

validateFcns(mpcobj_temp, rand(nx,1), 50);
clear mpcobj_temp p_temp stateFcn_temp outputFcn_temp;

%% Initialize storage for all patients
allResults = cell(numPatients, 1);

%% Extract patient data for parallel loop (parfor requires simple indexing)
patientIDs = patientData.ID;
patientBase = patientData.Base;
patientKtr = patientData.ktr;
patientSlope = patientData.slope;
patientGamma = patientData.gamma;

%% Run simulation for each patient in parallel
fprintf('\nStarting parallel simulation of %d patients...\n', numPatients);
tic; % Start timing

parfor patientIdx = 1:numPatients
    if ~isempty(getCurrentTask())
        fprintf('Processing Patient %s (Worker %d)\n', ...
            string(patientIDs(patientIdx)), getCurrentTask().ID);
    end
    
    % Set patient-specific parameters from CSV
    p = p_fixed;
    p.Base = patientBase(patientIdx);
    p.k_tr = patientKtr(patientIdx);
    p.slope = patientSlope(patientIdx);
    p.gamma = patientGamma(patientIdx);
    
    % Discrete-time plant for MPC (sample = 7 days)
    stateFcn = @(z,u) log( friberg_integrator(exp(z),u,p,Ts_wk) );
    outputFcn = @(x,u) exp(x(nx)); % measurement is ANC only
    
    % Build NL-MPC object
    mpcobj = nlmpc(nx,ny,nu);
    mpcobj.Ts = Ts_wk;
    mpcobj.Model.StateFcn = stateFcn;
    mpcobj.Model.OutputFcn = outputFcn;
    mpcobj.Model.IsContinuousTime = false;
    mpcobj.PredictionHorizon = mpcHorizon;
    mpcobj.ControlHorizon = mpcHorizon;
    
    % Constraints
    mpcobj.MV.Min = 0; % mg/m^2/d
    mpcobj.MV.Max = mpcMaxAllowableDose;
    
    % Cost function weights
    mpcobj.Weights.OutputVariables = mpcAncWeight;
    mpcobj.Weights.ManipulatedVariables = mpcDoseWeight;
    mpcobj.Weights.ManipulatedVariablesRate = 0;
    mpcobj.Weights.ECR = mpcECR;
    
    % Build Unscented Kalman Filter (state estimation only)
    ukf = unscentedKalmanFilter( ...
        stateFcn, ...
        @(x) exp(x(nx)), ...
        'HasAdditiveMeasurementNoise', true);
    
    ukf.Alpha = 1e-2;
    ukf.Beta = 2;
    ukf.Kappa = 0;
    
    ukf.MeasurementNoise = ancMeasNoiseSD.^2;
    ukf.ProcessNoise = diag(ukfStateProcessNoiseCovs);
    
    % Initialization
    x0 = [0.1; p.Base*p.k_ma/p.k_tr*ones(4,1); p.Base];
    
    zhat0 = [0.1; p.Base*p.k_ma/p.k_tr*ones(4,1); p.Base];
    ukf.State = log(zhat0);
    ukf.StateCovariance = diag(ukfInitialStateSDs.^2);
    
    % Closed-loop simulation
    xtrue = x0;
    dose_prev = 0;
    
    hist = struct();
    hist.x = nan(Nk,nx);
    hist.zhat = nan(Nk,nx);
    hist.y = nan(Nk,1);
    hist.u = nan(Nk,1);
    hist.t = Ts_wk*(0:Nk-1);
    hist.patientID = patientIDs(patientIdx);
    hist.trueParams = p;
    
    % Set reproducible random seed based on patient index
    rng(patientIdx); % Different seed for each patient, but reproducible
    
    for k = 1:Nk
        hist.x(k,:) = xtrue;
        
        % Measurement (clinic visit)
        % apply measurement noise in the log domain:
        yk = exp( log(xtrue(nx)) + ancMeasNoiseSD*randn );
        hist.y(k) = yk;
        
        % UKF update with 7-day batch
        try
            correct(ukf, yk);
        catch ME
            disp("Error in correct: " + ME.message);
            disp("    Patient ID = " + hist.patientID);
            disp("    k = " + k);
            break
        end
        
        hist.zhat(k,:) = ukf.State;
        
        % NMPC optimisation
        if yk < ancSafeLowerLim
            uk = 0;
        else
            try
                [uk,info] = nlmpcmove(mpcobj, hist.zhat(k,:), dose_prev, mpcAncTarget);
            catch ME
                disp("Error in nlmpcmove: " + ME.message);
                disp("    Patient ID = " + hist.patientID);
                disp("    k = " + k);
                break
            end
        end
        
        % Store & apply weekly constant dose
        hist.u(k) = uk;
        
        % Simulate true dynamics over next week
        xtrue = friberg_integrator(xtrue,uk,p,Ts_wk)
        
        % UKF time-update for the whole week (Ts_wk) with same uk
        try
            predict(ukf, uk);
        catch ME
            disp("Error in predict: " + ME.message);
            disp("    Patient ID = " + hist.patientID);
            disp("    k = " + k);
            break
        end
        
        dose_prev = uk;
    end
    
    % Store results for this patient
    allResults{patientIdx} = hist;
end

simTime = toc; % End timing
fprintf('\nAll simulations completed in %.2f minutes!\n', simTime/60);
fprintf('Average time per patient: %.2f seconds\n', simTime/numPatients);

%% Performance analysis (compare with serial execution time estimate)
% Rough estimate of serial time (based on typical single patient simulation time)
estimatedSerialTime = numPatients * 180; % Assume ~3 mins per patient
actualSpeedup = estimatedSerialTime / simTime;
if poolobj.Connected
    theoreticalMaxSpeedup = poolobj.NumWorkers;
else
    theoreticalMaxSpeedup = NaN;
end

fprintf('Performance Analysis:\n');
fprintf('  Workers used: %d\n', poolobj.NumWorkers);
fprintf('  Estimated speedup: %.1fx\n', actualSpeedup);
fprintf('  Theoretical max speedup: %.1fx\n', theoreticalMaxSpeedup);
fprintf('  Parallel efficiency: %.1f%%\n', 100 * actualSpeedup / theoreticalMaxSpeedup);

%% Save results
save("results/joint_ukf_no_adaptation_" + string(datetime("now"),"yyyyMMdd_hhmm"))

%% (Optional) Load results from a previous run

clear;

files = dir("results/joint_ukf_no_adaptation_*.mat");
if ~isempty(files)
    [~, idx] = max([files.datenum]);
    latestFile = fullfile("results", files(idx).name);
    fprintf('Loading most recent results: %s\n', latestFile);
    load(latestFile);
else
    fprintf('No results files found in results directory\n');
end

%% Close parallel pool
% delete(poolobj);

%% Summary statistics
fprintf('\n=== SUMMARY STATISTICS ===\n');
fprintf('Number of patients: %d\n', numPatients);

validResults = allResults(cellfun(@isstruct, allResults));

% Final doses
finalDoses = cellfun(@(x) x.u(end), validResults);
fprintf('Final doses (mg/m^2/d): Mean=%.2f, Std=%.2f, Min=%.2f, Max=%.2f\n', ...
    mean(finalDoses,'omitmissing'), std(finalDoses,'omitmissing'), min(finalDoses), max(finalDoses));

% Final ANC values
finalANC = cellfun(@(x) x.x(end,nx), validResults);
fprintf('Final ANC (G/L): Mean=%.2f, Std=%.2f, Min=%.2f, Max=%.2f\n', ...
    mean(finalANC,'omitmissing'), std(finalANC,'omitmissing'), min(finalANC), max(finalANC));

% Max ANC values
maxANC = cellfun(@(x) max(x.x(:,nx)), validResults);
failedTreatmentIdx = find(maxANC > 5);
fprintf('Max ANC was > 5 for %d patients\n', length(failedTreatmentIdx));

% Patients reaching target (within 10% of target ANC)
targetTolerance = 0.1 * mpcAncTarget;
patientsOnTarget = sum(abs(finalANC - mpcAncTarget) <= targetTolerance);
fprintf('Patients within 10%% of target ANC: %d/%d (%.1f%%)\n', ...
    patientsOnTarget, numPatients, 100*patientsOnTarget/numPatients);

% Find failed simulations
expectedFinalTime = Ts_wk * (Nk - 1);
failedSimIdx = [];
for i = 1:numPatients
    if isempty(allResults{i})
        actualTime = 0;
    else
        actualTime = max(allResults{i}.t(~isnan(allResults{i}.x(:,1))));
    end
    if actualTime < expectedFinalTime
        failedSimIdx = [failedSimIdx; i];
    end
end
fprintf('Failed simulations: %d/%d\n', length(failedSimIdx), numPatients);

%% Plotting results for all patients

% Generate distinct colors for each patient
colors = lines(numPatients);

% Plot 1: Drug doses and ANC for all patients
figure('Position', [100 100 800 600]);
subplot(2,1,1);
hold on;
for i = 1:numPatients
    if isempty(allResults{i})
        continue
    end
    [xs,ys] = stairs(allResults{i}.t, allResults{i}.u);
    ys(end) = NaN;
    plot(xs,ys, 'LineWidth', 1, ...
      'Color', [colors(i,:) 0.4], ...
      'DisplayName', sprintf('Patient %s', string(allResults{i}.patientID)));
end
xlim([0 max(allResults{1}.t)])
xlabel('Day');
ylabel('6-MP dose (mg/m^2/d)');
title('Drug Doses for All Patients');
legend('Location', 'eastoutside');
grid on;

subplot(2,1,2);
hold on;
for i = 1:numPatients
    if isempty(allResults{i})
        continue
    end
    % Plot true values:
    [xs,ys] = stairs(allResults{i}.t, allResults{i}.x(:,nx));
    ys(end) = NaN;
    plot(xs,ys, 'LineWidth', 1, ...
      'Color', [colors(i,:) 0.4], ...
      'DisplayName', sprintf('Patient %s', string(allResults{i}.patientID)));
    
    % Plot measured values:
    % plot(allResults{i}.t, allResults{i}.y, 'LineWidth', 1, ...
    %   'Color', [colors(i,:) 0.4], 'LineStyle', 'None', 'Marker', 'o', ...
    %   'DisplayName', sprintf('Patient %s', string(allResults{i}.patientID)));
end
yline(mpcAncTarget, 'b--', 'Target', 'LineWidth', 2);
xlim([0 max(allResults{1}.t)])
% ylim([0 3]);
xlabel('Day');
ylabel('ANC (10^9 cells/L)');
title('ANC Response for All Patients');
legend('Location', 'eastoutside');
grid on;

%% Plot true states in log domain

figure('Position',[100 100 600 900]);
tiledlayout(nx,1,'TileSpacing','tight');

for j = 1:nx
    nexttile; hold on; grid on;
    
    for i = 1:numPatients
        if isempty(allResults{i}), continue, end
        
        t = allResults{i}.t; % time vector
        plot(t, log(allResults{i}.x(:,j)), ...
            'LineWidth',1, ...
            'Color',[colors(i,:) 0.4], ...
            'DisplayName',sprintf('Patient %s',string(allResults{i}.patientID)));
    end

    xlim([0, max(cellfun(@(r)max(r.t), allResults(~cellfun(@isempty,allResults))))]);
    % ylim([-110 500])
    ylabel(sprintf('log(x_%d)', j));
    if j == nx
        xlabel('Time (days)')
    end
end

sgtitle('State Trajectories', 'FontSize', 14, 'FontWeight', 'bold');

%% Calculate relative state estimation errors

stateEstErrors = cell(numPatients,1);
colors = lines(numPatients);

for i = 1:numPatients
    if isempty(allResults{i})
        fprintf('Skipping patient %d (no results)\n', i);
        continue
    end
    
    % True and estimated states
    trueStates = allResults{i}.x;
    estimatedStates = exp(allResults{i}.zhat);
    
    % Calculate relative error
    stateEstErrors{i} = (estimatedStates - trueStates) ./ trueStates * 100;
end

%% Plot relative state estimation errors

figure('Position',[100 100 600 900]);
tiledlayout(nx,1,'TileSpacing','tight');

for j = 1:nx
    nexttile; hold on; grid on;
    
    for i = 1:numPatients
        if isempty(allResults{i}), continue, end
        
        t = allResults{i}.t; % time vector  
        plot(t, stateEstErrors{i}(:,j), ...
            'LineWidth',1, ...
            'Color',[colors(i,:) 0.4], ...
            'DisplayName',sprintf('Patient %s',string(allResults{i}.patientID)));
    end

    xlim([0, max(cellfun(@(r)max(r.t), allResults(~cellfun(@isempty,allResults))))]);
    ylim([-110 500])
    ylabel(sprintf('%% error in x_%d', j));
    if j == nx
        xlabel('Time (days)')
    end
end

sgtitle('State Estimation Errors Over Time', 'FontSize', 14, 'FontWeight', 'bold');

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
        x = max(x,1e-8); % Ensure positive states after update
    end
    xnext = x;
end

% PK/PK dynamics
function dx = friberg_rhs(x,dose,p)
    x = max(x,1e-8); % Ensure positive states

    C=x(1); P=x(2); T1=x(3); T2=x(4); T3=x(5); N=x(6);

    Edrug = p.slope*C;
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
function xnew = rk4(f,x,h)
    k1 = f(x);
    k2 = f(x + 0.5*h*k1);
    k3 = f(x + 0.5*h*k2);
    k4 = f(x + h*k3);
    xnew = x + h/6*(k1 + 2*k2 + 2*k3 + k4);
end

