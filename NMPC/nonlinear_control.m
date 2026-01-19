%% Nonlinear model predictive control using Dual Unscented Kalman Filter
% Simplified pharmacokinetics compared to Jost 2020
% MODIFIED: Replaced Joint UKF with Dual UKF for estimation.

clear; clc; close all;

userdir = string(java.lang.System.getProperty('user.home')) + "/";
codeFolder = userdir + "UofT/Research/Leukemia/Code/";

% Path to CSV file with patient parameters (Data Sheet 1 from Jost 2020 supplementary materials)
csvFilePath = codeFolder + "Jost 2020 Data Analysis/Supplementary Materials/Data Sheet 1.CSV";

%% Read patient parameters from CSV
patientData = readtable(csvFilePath);
numPatients = height(patientData);

fprintf('Found %d patients in CSV file\n', numPatients);

%% Parallel pool setup
% Check if parallel pool is already running
% poolobj = gcp('nocreate');
% if isempty(poolobj)
%     fprintf('Starting parallel pool...\n');
%     parpool('local'); % Start local parallel pool
%     poolobj = gcp;
% end
% fprintf('Using parallel pool with %d workers\n', poolobj.NumWorkers);

%% Fixed model/PK parameters (same for all patients)
p_fixed.ke = 0.22; % 1/day elimination of 6-MP surrogate concentration
p_fixed.F_over_V = 0.00033; % lumped oral bioavailability / volume term
p_fixed.k_ma = 2.3765; % 1/day

% Parameter estimates (only Base, k_tr, slope, and gamma are adapted)
phat_fixed.ke = p_fixed.ke;
phat_fixed.F_over_V = p_fixed.F_over_V;
phat_fixed.k_ma = 2.3765; % 1/day

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
mpcMaxAllowableDose = 200; % Maximum 6-MP dose (mg/m^2)
mpcAncWeight = 1; % Weight for ANC
mpcDoseWeight = 0; % Weight for 6-MP dose
mpcHorizon = 8; % Number of weeks to predict and control
mpcECR = 1e4; % Slack variable tuning weight

% UKF parameters (Important: these are in the log domain)
ukfStateProcessNoiseCovs = [1e-6 1e-4 1e-4 1e-4 1e-4 1e-6];
ukfParamProcessNoiseCovs = [1e-5 1e-5 1e-5 1e-4];
ukfInitialStateSDs = [1e-3 0.125304565*ones(1,4) 0.103687761];
ukfInitialParamSDs = [0.103687761 0.059575249 0.161059312 0.044632609]; % SD of log-params from Jost Data Sheet 1

%% Initialize storage for all patients
allResults = cell(numPatients, 1);

%% Extract patient data for parallel loop (parfor requires simple indexing)
patientIDs = patientData.ID;
patientBase = patientData.Base;
patientKtr = patientData.ktr;
patientSlope = patientData.slope;
patientGamma = patientData.gamma;

% Extra margin (as fraction) to allow for possible parameter estimates
% e.g., paramBoundMargin=0.1 means 10% above highest and below lowest
% population parameter values
paramBoundMargin = 0.1;
paramSafeLowerLim = [min(patientBase); min(patientKtr); min(patientSlope); min(patientGamma)] .* (1-paramBoundMargin);
paramSafeUpperLim = [max(patientBase); max(patientKtr); max(patientSlope); max(patientGamma)] .* (1+paramBoundMargin);

% Leave-one-out averages
sum_Base = sum(patientBase);
sum_k_tr = sum(patientKtr);
sum_slope = sum(patientSlope);
sum_gamma = sum(patientGamma);

avg_Base_excl = (sum_Base - patientBase) ./ (numPatients - 1);
avg_ktr_excl = (sum_k_tr - patientKtr) ./ (numPatients - 1);
avg_slope_excl = (sum_slope - patientSlope) ./ (numPatients - 1);
avg_gamma_excl = (sum_gamma - patientGamma) ./ (numPatients - 1);

%% Run simulation for each patient in parallel
fprintf('\nStarting parallel simulation of %d patients...\n', numPatients);
tic; % Start timing

for patientIdx = 2%:numPatients
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
    
    % Set patient-specific parameter *estimates*
    phat = phat_fixed;
    phat.Base = avg_Base_excl(patientIdx);
    phat.k_tr = avg_ktr_excl(patientIdx);
    phat.slope = avg_slope_excl(patientIdx);
    phat.gamma = avg_gamma_excl(patientIdx);
    
    % Discrete-time plant for MPC (sample = 7 days)
    % This model is for the MPC's internal prediction, it uses the augmented state.
    mpcStateFcn = @(z,u) [
        log( friberg_integrator(exp(z(1:nx)),u,makeParamStruct(p,exp(z(nx+1:nx+np))),Ts_wk) );
        z(nx+1:end) % parameters assumed are constant for MPC prediction horizon
    ];
    mpcOutputFcn = @(x,u) exp(x(nx)); % measurement is ANC only
    
    % Build NL-MPC object
    mpcobj = nlmpc(nx+np,ny,nu);
    mpcobj.Ts = Ts_wk;
    mpcobj.Model.StateFcn = mpcStateFcn;
    mpcobj.Model.OutputFcn = mpcOutputFcn;
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
    
    % ===== DUAL UKF SETUP =====
    
    % 1. State UKF
    stateUkf = unscentedKalmanFilter(...
        @(x_log, u, p_struct) state_transition_dukf(x_log, u, p_struct, Ts_wk), ...
        @(x_log) state_measurement_dukf(x_log, nx), ...
        'HasAdditiveMeasurementNoise', true);
    stateUkf.Alpha = 1e-2;
    stateUkf.Beta = 2;
    stateUkf.Kappa = 0;
    stateUkf.MeasurementNoise = ancMeasNoiseSD.^2;
    stateUkf.ProcessNoise = diag(ukfStateProcessNoiseCovs);

    % 2. Parameter UKF
    paramUkf = unscentedKalmanFilter(...
        @(p_log) p_log, ... % Parameter state transition is identity
        @(p_log, x_log, u, p_fix) param_measurement_dukf(p_log, x_log, u, p_fix, Ts_wk, nx), ...
        'HasAdditiveMeasurementNoise', true);
    paramUkf.Alpha = 1e-2;
    paramUkf.Beta = 2;
    paramUkf.Kappa = 0;
    % The measurement noise for the parameter filter reflects the uncertainty in the measurement prediction
    paramUkf.MeasurementNoise = ancMeasNoiseSD.^2;
    paramUkf.ProcessNoise = diag(ukfParamProcessNoiseCovs);
    
    % Initialization
    x0 = [0.1; p.Base*p.k_ma/p.k_tr*ones(4,1); p.Base];
    
    xhat0 = [1e-3; phat.Base*phat.k_ma/phat.k_tr*ones(4,1); phat.Base];
    phat0 = [phat.Base; phat.k_tr; phat.slope; phat.gamma];
    
    stateUkf.State = log(xhat0);
    stateUkf.StateCovariance = diag(ukfInitialStateSDs.^2);
    
    paramUkf.State = log(phat0);
    paramUkf.StateCovariance = diag(ukfInitialParamSDs.^2);
    
    % Closed-loop simulation
    xtrue = x0;
    dose_prev = 0;
    
    hist = struct();
    hist.x = nan(Nk,nx);
    hist.zhat = nan(Nk,nx+np); % Storing log([x_hat; p_hat])
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
        
        % ===== DUAL UKF Update Step =====
        
        % Get prior estimates needed for corrections
        xhat_log_prior = stateUkf.State;
        
        try
            % 1. Correct State UKF
            % The measurement function for the state filter is h(x)=exp(x(nx))
            % and does not require additional parameters.
            correct(stateUkf, yk);
            
            % 2. Correct Parameter UKF
            % The measurement function for the parameter filter needs the prior
            % state estimate and the previous input to predict the measurement yk.
            correct(paramUkf, yk, xhat_log_prior, dose_prev, p_fixed);

        catch ME
            disp("Error in UKF correct step: " + ME.message);
            disp("    Patient ID = " + hist.patientID);
            disp("    k = " + k);
            break
        end

        % Apply safety constraints on parameters
        currentParamEst = exp(paramUkf.State);
        currentParamEst(currentParamEst > paramSafeUpperLim) = paramSafeUpperLim(currentParamEst > paramSafeUpperLim);
        currentParamEst(currentParamEst < paramSafeLowerLim) = paramSafeLowerLim(currentParamEst < paramSafeLowerLim);
        try
            paramUkf.State = log(currentParamEst);
        catch ME
            disp("Error when applying parameter constraints: " + ME.message);
            disp("    Patient ID = " + hist.patientID);
            disp("    k = " + k);
            break
        end
        
        % Store combined log-state for history and MPC
        zhat_log_k = [stateUkf.State; paramUkf.State];
        hist.zhat(k,:) = zhat_log_k;
        
        % NMPC optimisation
        if yk < ancSafeLowerLim
            uk = 0;
        else
            try
                [uk,info] = nlmpcmove(mpcobj, zhat_log_k, dose_prev, mpcAncTarget);
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
        xtrue = friberg_integrator(xtrue,uk,p,Ts_wk);
        
        % ===== DUAL UKF Predict Step =====
        
        % Get latest parameter estimate to use in state prediction
        phat_posterior_log = paramUkf.State;
        phat_posterior_struct = makeParamStruct(phat_fixed, exp(phat_posterior_log));
        
        try
            % 1. Predict State UKF
            % The state transition f(x, u, p) needs the latest parameter estimate.
            predict(stateUkf, uk, phat_posterior_struct);
            
            % 2. Predict Parameter UKF
            % The parameter transition is f(p) = p, which needs no extra inputs.
            predict(paramUkf);
        catch ME
            disp("Error in UKF predict step: " + ME.message);
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

%% Save results
save("results/dual_ukf_" + string(datetime("now"),"yyyyMMdd_hhmm"))

%% (Optional) Load results from a previous run

clear;

files = dir("results/dual_ukf_*.mat");
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
% xlim([0 max(allResults{1}.t)])
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
    plot(allResults{i}.t, allResults{i}.y, 'LineWidth', 1, ...
      'Color', [colors(i,:) 0.4], 'LineStyle', 'None', 'Marker', 'o', ...
      'DisplayName', sprintf('Patient %s', string(allResults{i}.patientID)));
end
yline(mpcAncTarget, 'b--', 'Target', 'LineWidth', 2);
% xlim([0 max(allResults{1}.t)])
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
             'LineWidth', 1, 'Color', [colors(i,:) 0.6], ...
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

% === DUAL UKF SPECIFIC FUNCTIONS ===

% DUAL UKF - State Transition Function
function xnext_log = state_transition_dukf(x_log, u, p_struct, Ts)
    % Unpack states from log domain
    x = exp(x_log);
    % Integrate dynamics
    xnext = friberg_integrator(x, u, p_struct, Ts);
    % Repack states into log domain
    xnext_log = log(xnext);
end

% DUAL UKF - State Measurement Function
function y = state_measurement_dukf(x_log, nx)
    % Predicts measurement (ANC) from the log-state vector
    y = exp(x_log(nx));
end

% DUAL UKF - Parameter "Measurement" Function
function y_pred = param_measurement_dukf(p_log, x_log_prior, u, p_fixed, Ts, nx)
    % This function predicts the measurement 'y' based on a given set of
    % parameters (from sigma points), the prior state estimate, and the input.

    % Unpack parameters from log domain and create struct
    p_vec = exp(p_log);
    p_struct = makeParamStruct(p_fixed, p_vec);
    
    % Unpack prior state estimate from log domain
    x_prior = exp(x_log_prior);
    
    % Predict the next state using these parameters
    x_pred = friberg_integrator(x_prior, u, p_struct, Ts);
    
    % The predicted measurement is the ANC from the predicted state
    % Output is in linear space, as required by the filter's residual calculation.
    y_pred = x_pred(nx);
end


% === GENERAL MODEL AND UTILITY FUNCTIONS ===

% Integrated Friberg dynamics over the time period Delta
function xnext = friberg_integrator(x,u,p,Delta)
    % Integrate PK+PD over Delta (days) with daily bolus dose u.
    dt = 1/24; % Internal step size
    Nt = round(Delta/dt);
    for n = 1:Nt
        drugRate = 0;
        % Dose event happens once every day at hour zero
        if mod((n-1)*dt,1) == 0
            drugRate = u/dt; % Normalized by dt so total drug input is u
        end
        x = rk4(@(z) friberg_rhs(z,drugRate,p), x, dt);
        x = max(x,1e-8); % Ensure positive states after update
    end
    xnext = x;
end

% PK/PK dynamics
function dx = friberg_rhs(x,drugRate,p)
    C=x(1); P=x(2); T1=x(3); T2=x(4); T3=x(5); N=x(6);

    Edrug = p.slope*C;
    feedback = (p.Base/N)^p.gamma;

    % Single compartment PK:
    dC  = -p.ke*C + p.F_over_V*drugRate;

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

function pstruct = makeParamStruct(p,pvec)
    pstruct.ke = p.ke;
    pstruct.F_over_V = p.F_over_V;
    pstruct.Base = pvec(1);
    pstruct.k_tr = pvec(2);
    pstruct.slope = pvec(3);
    pstruct.gamma = pvec(4);
    pstruct.k_ma = p.k_ma;
end
