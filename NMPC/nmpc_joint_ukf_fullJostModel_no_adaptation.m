%% Nonlinear model predictive control using Unscented Kalman Filter
% Model from Jost 2020
% No parameter adaptation. Parameters are assumed fixed and known.

clear; clc; close all;

% Path to CSV file with patient parameters (Data Sheet 1 from Jost 2020 supplementary materials)
paramFilePath = "C:\Users\trent\UofT\Research\Leukemia\Code\Jost 2020 Data Analysis\Supplementary Materials\Data Sheet 1.CSV";
dataFilePath = "C:\Users\trent\UofT\Research\Leukemia\Code\Jost 2020 Data Analysis\Supplementary Materials\Data Sheet 2.CSV";

%% Read patient parameters from CSV
patientData = readtable(paramFilePath);
numPatients = height(patientData);

fprintf('Found %d patients in CSV file\n', numPatients);

%% Read patient treatment data from CSV
T = readtable(dataFilePath);

% Extract unique BSA values for each patient
patientBSA = zeros(numPatients, 1);
for i = 1:numPatients
    patientID = patientData.ID(i);
    % Find rows in treatment data that match this patient ID
    matchingRows = T.ID == patientID;
    if any(matchingRows)
        % BSA should be the same for all rows of the same patient
        patientBSA(i) = T.BSA(find(matchingRows, 1));
    else
        warning('No BSA data found for patient %s, using default value', string(patientID));
        patientBSA(i) = 1.0; % Default BSA value
    end
end

fprintf('Extracted BSA values for %d patients\n', numPatients);

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
p_fixed.F = 0.22; % Oral bioavailability
p_fixed.k_a = 31.2; % Absorption rate constant of 6MP
p_fixed.k_20 = 12.72; % Elimination rate constant of 6MP
p_fixed.FM_3 = 0.019; % Fractional metabolic transformation into 6-TGN
p_fixed.k_me = 9.9216; % Metabolic transformation rate constant of 6MP -> 6-GTN or 6-MPN
p_fixed.CL_6tgn = 0.219*1^1.16; % BSA-dependent clearance of 6-TGN (placeholder)
p_fixed.k_ma = 2.3765; % Death rate of matured neutrophils

% Parameter estimates (same as p_fixed)
phat_fixed.F = p_fixed.F;
phat_fixed.k_a = p_fixed.k_a;
phat_fixed.k_20 = p_fixed.k_20;
phat_fixed.FM_3 = p_fixed.FM_3;
phat_fixed.k_me = p_fixed.k_me;
phat_fixed.CL_6tgn = p_fixed.CL_6tgn;
phat_fixed.k_ma = p_fixed.k_ma;

ancMeasNoiseSD = 0.1;

%% Simulation parameters
Ts_wk = 7; % days
nx = 8;
ny = 1;
nu = 1;
np = 4;
Nk = 30; % weeks to simulate

ancSafeLowerLim = 0.5; % If ANC is below this level, give zero dose

% MPC parameters
mpcAncTarget = 1.5; % target ANC
mpcMaxAllowableDose = 200; % Maximum 6-MP dose (mg/day)
mpcAncWeight = 1; % Weight for ANC
mpcDoseWeight = 0.01; % Weight for 6-MP dose
mpcHorizon = 8; % Number of weeks to predict and control
mpcECR = 1e5; % Slack variable tuning weight

% UKF parameters (Important: these are in the log domain)
ukfStateProcessNoiseCovs = [1e-6 1e-6 1e-6 1e-4 1e-4 1e-4 1e-4 1e-6];
ukfParamProcessNoiseCovs = [1e-5 1e-5 1e-5 1e-4];
ukfInitialStateSDs = [1e-2*ones(1,3) 0.125304565*ones(1,4) 0.103687761];
ukfInitialParamSDs = [0.103687761 0.059575249 0.161059312 0.044632609]; % SD of log-params from Jost Data Sheet 1

%% Validate functions once before parallel execution
fprintf('Validating MPC functions...\n');
p_temp = p_fixed;
p_temp.Base = 2.5;
p_temp.k_tr = 0.15;
p_temp.slope = 0.27;
p_temp.gamma = 0.8;

stateFcn_temp = @(z,u) log( jost_integrator(exp(z(1:nx)),u,p_temp,Ts_wk) );
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
% patientBSA is already extracted above

%% Run simulation for each patient in parallel
fprintf('\nStarting parallel simulation of %d patients...\n', numPatients);
tic; % Start timing

% patientIdSubset = [184,162,194,176,121,110];

parfor patientIdx = 82:numPatients
    % patientIdx = find(patientIdSubset(i) == patientIDs, 1)
    
    if ~isempty(getCurrentTask())
        fprintf('Processing Patient %s (Worker %d)\n', ...
            string(patientIDs(patientIdx)), getCurrentTask().ID);
    end
    
    % Set patient-specific *true* parameters from CSV
    p = p_fixed;
    p.Base = patientBase(patientIdx);
    p.k_tr = patientKtr(patientIdx);
    p.slope = patientSlope(patientIdx);
    p.gamma = patientGamma(patientIdx);
    p.CL_6tgn = 0.219 * patientBSA(patientIdx)^1.16;
    
    % Discrete-time plant for MPC (sample = 7 days)
    stateFcn = @(z,u) log( jost_integrator(exp(z(1:nx)),u,p,Ts_wk) );
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
    mpcobj.MV.Min = 0; % mg/day
    mpcobj.MV.Max = mpcMaxAllowableDose;
    
    % Cost function weights
    mpcobj.Weights.OutputVariables = mpcAncWeight;
    mpcobj.Weights.ManipulatedVariables = mpcDoseWeight;
    mpcobj.Weights.ManipulatedVariablesRate = 0;
    mpcobj.Weights.ECR = mpcECR;
    
    % Build Unscented Kalman Filter
    jointUkf = unscentedKalmanFilter( ...
        stateFcn, ...
        @(x) exp(x(nx)), ...
        'HasAdditiveMeasurementNoise', true);
    
    jointUkf.Alpha = 1e-2;
    jointUkf.Beta = 2;
    jointUkf.Kappa = 0;
    
    jointUkf.MeasurementNoise = ancMeasNoiseSD.^2;
    jointUkf.ProcessNoise = diag(ukfStateProcessNoiseCovs);
    
    % Initialization
    x0 = [1e-3*ones(3,1); p.Base*p.k_ma/p.k_tr*ones(4,1); p.Base]; % 8 elements
    
    zhat0 = [1e-3*ones(3,1); p.Base*p.k_ma/p.k_tr*ones(4,1); p.Base];
    jointUkf.State = log(zhat0);
    jointUkf.StateCovariance = diag(ukfInitialStateSDs.^2);
    
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
            correct(jointUkf, yk);
        catch ME
            disp("Error in correct: " + ME.message);
            disp("    Patient ID = " + hist.patientID);
            disp("    k = " + k);
            break
        end

        hist.zhat(k,:) = jointUkf.State;
        
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
        xtrue = jost_integrator(xtrue,uk,p,Ts_wk);
        
        % UKF time-update for the whole week (Ts_wk) with same uk
        try
            predict(jointUkf, uk);
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
save("results/joint_ukf_fullJostModel" + string(datetime("now"),"yyyyMMdd_hhmm"))

%% (Optional) Load results from a previous run

% clear; load("results/joint_ukf_20250611_0149.mat")
% clear; load("results/joint_ukf_20250619_1148.mat")
% clear; load("results/joint_ukf_20250619_0157.mat")
% clear; load("results/joint_ukf_20250620_0307.mat")
% clear; load("results/joint_ukf_20250626_0830.mat")
% clear; load("results/joint_ukf_fullJostModel20250628_0834.mat")
% clear; load("results/joint_ukf_fullJostModel20250715_0428.mat")
clear; load("results/joint_ukf_fullJostModel20250717_1206.mat")

%% Close parallel pool
% delete(poolobj);

%% Summary statistics
fprintf('\n=== SUMMARY STATISTICS ===\n');
fprintf('Number of patients: %d\n', numPatients);

validResults = allResults(cellfun(@isstruct, allResults));

% Final doses
finalDoses = cellfun(@(x) x.u(end), validResults);
fprintf('Final doses (mg/day): Mean=%.2f, Std=%.2f, Min=%.2f, Max=%.2f\n', ...
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
      'DisplayName', sprintf('Patient %s', string(allResults{i}.patientID)));
end
xlim([0 max(allResults{1}.t)])
xlabel('Day');
ylabel('6-MP dose (mg/day)');
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
    estimatedStates = exp(allResults{i}.zhat(:, 1:nx));
    
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

%% Plot parameter estimation errors (absolute)

paramNames = {'Base', 'k_{tr}', 'slope', '\gamma'};
errorFields = {'Base', 'k_tr', 'slope', 'gamma'};

figure('Position', [200 200 1000 800]);
for paramIdx = 1:np
    subplot(2,2,paramIdx);
    hold on;
    for i = 1:numPatients
        if isempty(allResults{i})
            continue
        end
        lineData = stairs(allResults{i}.t, exp(allResults{i}.zhat(:,nx+paramIdx)), ...
               'LineWidth', 1.2, ...
               'DisplayName', sprintf('Patient %s', string(allResults{i}.patientID)));
        % Add true value line for this patient
        trueValue = allResults{i}.trueParams.(errorFields{paramIdx});
        yline(trueValue, '--', 'LineWidth', 1, ...
              'HandleVisibility', 'off');
    end
    xlim([0 max(allResults{1}.t)])
    xlabel('Day');
    ylabel(paramNames{paramIdx});
    title(sprintf('Parameter Estimation: %s', paramNames{paramIdx}));
    legend('Location', 'eastoutside');
    grid on;
end

%% Calculate relative parameter estimation errors

paramNames = {'Base', 'k_{tr}', 'slope', '\gamma'};
errorFields = {'Base', 'k_tr', 'slope', 'gamma'};

for i = 1:numPatients
    if isempty(allResults{i})
        fprintf('Skipping patient %d (no results)\n', i);
        continue
    end
    
    % Get true parameter values for this patient
    trueParams = allResults{i}.trueParams;
    % Get estimated parameter trajectories (convert from log space)
    estimatedParams = exp(allResults{i}.zhat(:, nx+1:nx+np));
    
    % Calculate relative errors for each time step and parameter

    paramEstErrors.Base(i, :) = (estimatedParams(:, 1) - trueParams.Base) / trueParams.Base * 100;
    finalErrors.Base(i) = paramEstErrors.Base(i, end);

    paramEstErrors.k_tr(i, :) = (estimatedParams(:, 2) - trueParams.k_tr) / trueParams.k_tr * 100;
    finalErrors.k_tr(i) = paramEstErrors.k_tr(i, end);

    paramEstErrors.slope(i, :) = (estimatedParams(:, 3) - trueParams.slope) / trueParams.slope * 100;
    finalErrors.slope(i) = paramEstErrors.slope(i, end);
    
    % gamma parameter (index 4)
    paramEstErrors.gamma(i, :) = (estimatedParams(:, 4) - trueParams.gamma) / trueParams.gamma * 100;
    finalErrors.gamma(i) = paramEstErrors.gamma(i, end);
    
    % fprintf('Patient %s: Final relative errors - Base: %.1f%%, k_tr: %.1f%%, slope: %.1f%%, gamma: %.1f%%\n', ...
    %     string(allResults{i}.patientID), finalErrors.Base(i), finalErrors.k_tr(i), ...
    %     finalErrors.slope(i), finalErrors.gamma(i));
end

%% Plot relative parameter estimation errors over time

figure('Position', [300 300 1200 800]);
for paramIdx = 1:np
    subplot(2,2,paramIdx);
    hold on;
    
    % Plot error trajectories for all patients
    for i = 1:numPatients
        if isempty(allResults{i})
            continue
        end
        plot(allResults{i}.t, paramEstErrors.(errorFields{paramIdx})(i, :), ...
             'LineWidth', 1, ...
             'DisplayName', sprintf('Patient %s', string(allResults{i}.patientID)));
    end
    
    % Add zero error reference line
    yline(0, 'k--', 'Perfect Estimation', 'LineWidth', 2);
    
    xlim([0 max(allResults{1}.t)])
    xlabel('Day');
    ylabel('Relative Error (%)');
    title(sprintf('Parameter Estimation Error: %s', paramNames{paramIdx}));
    grid on;
    
    % Set reasonable y-axis limits
    allErrors = paramEstErrors.(errorFields{paramIdx});
    validErrors = allErrors(~isnan(allErrors));
    if ~isempty(validErrors)
        errorRange = [min(validErrors), max(validErrors)];
        yMargin = 0.1 * (errorRange(2) - errorRange(1));
        ylim([errorRange(1) - yMargin, errorRange(2) + yMargin]);
    end
end

sgtitle('Parameter Estimation Errors Over Time', 'FontSize', 14, 'FontWeight', 'bold');

%% Plot final parameter estimation error distributions

figure('Position', [400 400 1000 600]);
for paramIdx = 1:np
    subplot(2,2,paramIdx);
    
    validErrors = finalErrors.(errorFields{paramIdx})(~isnan(finalErrors.(errorFields{paramIdx})));
    
    if ~isempty(validErrors)
        histogram(validErrors, 'Normalization', 'probability', ...
                 'FaceColor', [0.7 0.7 0.9], 'EdgeColor', 'black', ...
                 'FaceAlpha', 0.7);
        
        % Add vertical line for mean error
        xline(mean(validErrors), 'r--', sprintf('Mean: %.1f%%', mean(validErrors)), ...
              'LineWidth', 2, 'LabelHorizontalAlignment', 'center');
        
        % Add vertical line for zero error
        xline(0, 'k--', 'Perfect', 'LineWidth', 2, ...
              'LabelHorizontalAlignment', 'center');
        
        xlabel('Final Relative Error (%)');
        ylabel('Probability');
        title(sprintf('Final %s Estimation Error Distribution', paramNames{paramIdx}));
        grid on;
    else
        text(0.5, 0.5, 'No valid data', 'HorizontalAlignment', 'center', ...
             'VerticalAlignment', 'middle', 'FontSize', 12);
        title(sprintf('Final %s Estimation Error Distribution', paramNames{paramIdx}));
    end
end

sgtitle('Distribution of Final Parameter Estimation Errors', 'FontSize', 14, 'FontWeight', 'bold');

%% Helper functions

% Integrated Jost 2020 dynamics over the time period Delta
function xnext = jost_integrator(x,u,p,Delta)
    % Integrate PK+PD over Delta (days) with daily bolus dose u.
    dt = 1/24; % internal step of 1 hour
    Nt = round(Delta/dt);
    for n = 1:Nt
        x = rk4(@(z) jost_rhs(z,u,p), x, dt);
        x = max(x,1e-8); % Ensure positive states after update
    end
    xnext = x;
end

% PK/PK dynamics
function dx = jost_rhs(x,dose,p)
    x = max(x,1e-8); % Ensure positive states

    Edrug = p.slope*x(3);
    feedback = (p.Base/x(8))^p.gamma;

    dx = zeros(8,1);

    % PK from Jost 2020
    dx(1) = -p.k_a * x(1) + p.F * dose;
    dx(2) = p.k_a * x(1) - p.k_20 * x(2);
    dx(3) = p.FM_3*p.k_me * x(2) - p.CL_6tgn * x(3);

    % Friberg model with 3 transitional states:
    dx(4) = p.k_tr*x(4)*(1-Edrug)*feedback - p.k_tr*x(4);
    dx(5) = p.k_tr*x(4) - p.k_tr*x(5);
    dx(6) = p.k_tr*x(5) - p.k_tr*x(6);
    dx(7) = p.k_tr*x(6) - p.k_tr*x(7);
    dx(8) = p.k_tr*x(7) - p.k_ma*x(8);
end

% 4th order Runge-Kutta
function xnew = rk4(f,x,h)
    k1 = f(x);
    k2 = f(x + 0.5*h*k1);
    k3 = f(x + 0.5*h*k2);
    k4 = f(x + h*k3);
    xnew = x + h/6*(k1 + 2*k2 + 2*k3 + k4);
end
