%% Robust Nonlinear model predictive control using Joint UKF and Genetic Algorithm
% Model from Jost 2020.
% Robustness achieved by optimizing expected cost over parameter uncertainty using Unscented Transform.

clear; clc; close all;

% --- Setup Paths ADJUST AS NEEDED ---
userdir = string(java.lang.System.getProperty('user.home')) + "/";
codeFolder = userdir + "UofT/Research/Leukemia/Code/";

% Define paths to the Jost 2020 Data Sheets 1 (parameters) and 2 (BSA data).
paramFilePath = codeFolder + "Jost 2020 Data Analysis/Supplementary Materials/Data Sheet 1.CSV";
dataFilePath = codeFolder + "Jost 2020 Data Analysis/Supplementary Materials/Data Sheet 2.CSV";

if ~exist('paramFilePath', 'var') || ~exist('dataFilePath', 'var')
    error('Please define paramFilePath and dataFilePath variables pointing to the data sheets.');
end

% Create results directory
if ~exist('results', 'dir')
   mkdir('results')
end

%% Read patient parameters and BSA data
try
    patientData = readtable(paramFilePath);
    numPatients = height(patientData);
    fprintf('Found %d patients in parameter file.\n', numPatients);
catch
    error('Could not read parameter file: %s', paramFilePath);
end

try
    T = readtable(dataFilePath);
    patientBSA = ones(numPatients, 1); % Default BSA
    for i = 1:numPatients
        patientID = patientData.ID(i);
        matchingRows = T.ID == patientID;
        if any(matchingRows)
            patientBSA(i) = T.BSA(find(matchingRows, 1));
        end
    end
    fprintf('Extracted BSA values.\n');
catch
     warning('Could not read treatment data file. Proceeding with default BSA=1.0. Path: %s', dataFilePath);
     patientBSA = ones(numPatients, 1);
end

%% Parallel pool setup
poolobj = gcp('nocreate');
if isempty(poolobj)
    fprintf('Starting parallel pool...\n');
    try
        parpool('local');
        poolobj = gcp;
    catch
        warning('Failed to start parallel pool. Proceeding sequentially.');
        poolobj = [];
    end
end

if ~isempty(poolobj)
    fprintf('Using parallel pool with %d workers\n', poolobj.NumWorkers);
    useParallel = true;
    parforArg = Inf;
else
    useParallel = false;
    parforArg = 0;
end

%% Fixed model/PK parameters
p_fixed.F = 0.22; p_fixed.k_a = 31.2; p_fixed.k_20 = 12.72;
p_fixed.FM_3 = 0.019; p_fixed.k_me = 9.9216;
p_fixed.CL_6tgn = nan; % Overwritten based on BSA
p_fixed.k_ma = 2.3765;

phat_fixed = p_fixed; % Structure for fixed estimates

%% Simulation and Controller parameters
Ts_wk = 7; % days (sampling time)
nx = 8; % Number of states
np = 4; % Number of estimated parameters [Base, k_tr, slope, gamma]
Nk = 30; % weeks to simulate

ancTarget = 1.5;
% ancSafeLowerLim = 0.2;
ancMeasNoiseSD = 0.1; % SD in log domain

maxAllowableDose = 200;
ancWeight = 1e3;
doseWeight = 5e-2;
predictionHorizon = 6;

possible_doses = 0:25:maxAllowableDose;
n_possible_doses = length(possible_doses);

% UKF Noise Covariances (Log domain)
ukfStateProcessNoiseCovs = [1e-6 1e-6 1e-6 5e-4 5e-4 5e-4 5e-4 1e-6];
ukfParamProcessNoiseCovs = [1e-4 1e-4 1e-4 1e-4];
ukfInitialStateSDs = [1e-2*ones(1,3) 0.125304565*ones(1,4) 0.103687761];
ukfInitialParamSDs = [0.103687761 0.059575249 0.161059312 0.044632609];

% Unscented Transform Parameters (for UKF estimation and Robust MPC calculation)
ut_params.Alpha = 1e-2; % Controls the spread of the sigma points
ut_params.Beta = 2;     % Optimal for Gaussian distributions
ut_params.Kappa = 0;    % Secondary scaling parameter

%% Initialization and Data Extraction
allResults = cell(numPatients, 1);

patientIDs = patientData.ID;
patientBase = patientData.Base;
patientKtr = patientData.ktr;
patientSlope = patientData.slope;
patientGamma = patientData.gamma;

% Calculate Leave-one-out averages for initial estimates
sum_log_Base = sum(log(patientBase));
sum_log_k_tr = sum(log(patientKtr));
sum_log_slope = sum(log(patientSlope));
sum_log_gamma = sum(log(patientGamma));

avg_Base_excl = exp( (sum_log_Base - log(patientBase)) ./ (numPatients - 1) );
avg_ktr_excl = exp( (sum_log_k_tr - log(patientKtr)) ./ (numPatients - 1) );
avg_slope_excl = exp( (sum_log_slope - log(patientSlope)) ./ (numPatients - 1) );
avg_gamma_excl = exp( (sum_log_gamma - log(patientGamma)) ./ (numPatients - 1) );

%% Run simulation for each patient
fprintf('\nStarting simulation of %d patients using Robust MPC (UT, Parameter Uncertainty)...\n', numPatients);
tic;

% parfor (patientIdx = 1:numPatients, parforArg)
for patientIdx = 1:numPatients % Use for debugging
    
    if useParallel && ~isempty(getCurrentTask())
        fprintf('Processing Patient %s (Worker %d)\n', string(patientIDs(patientIdx)), getCurrentTask().ID);
    else
        fprintf('Processing Patient %s\n', string(patientIDs(patientIdx)));
    end
    
    % 1. Setup Patient Parameters (True and Initial Estimates)
    p = p_fixed;
    p.Base = patientBase(patientIdx);
    p.k_tr = patientKtr(patientIdx);
    p.slope = patientSlope(patientIdx);
    p.gamma = patientGamma(patientIdx);
    p.CL_6tgn = 0.219 * patientBSA(patientIdx)^1.16;
    
    phat = phat_fixed;
    % Initial estimates for varying parameters
    phat.Base = avg_Base_excl(patientIdx);
    phat.k_tr = avg_ktr_excl(patientIdx);
    phat.slope = avg_slope_excl(patientIdx);
    phat.gamma = avg_gamma_excl(patientIdx);
    phat.CL_6tgn = p.CL_6tgn; % Considered known
    
    % 2. Setup UKF
    % State transition function (required by MATLAB UKF object)
    % 'phat' captures the fixed parameters (like CL_6tgn) for this patient.
    stateFcn = @(z,u) [
        log( jost_integrator(exp(z(1:nx)),u,makeParamStruct(phat,exp(z(nx+1:nx+np))),Ts_wk) );
        z(nx+1:end) % Parameters assumed constant during short prediction step
    ];
    
    jointUkf = unscentedKalmanFilter(stateFcn, @(z) z(nx), 'HasAdditiveMeasurementNoise', true);
    
    % Set UT parameters for the filter
    jointUkf.Alpha = ut_params.Alpha;
    jointUkf.Beta = ut_params.Beta;
    jointUkf.Kappa = ut_params.Kappa;
    
    jointUkf.MeasurementNoise = ancMeasNoiseSD.^2;
    jointUkf.ProcessNoise = blkdiag( diag(ukfStateProcessNoiseCovs), diag(ukfParamProcessNoiseCovs) );
    
    % Initialization
    % True initial state
    x0 = [1e-3*ones(3,1); p.Base*p.k_ma/p.k_tr*ones(4,1); p.Base];
    
    % Estimated initial state and parameters
    zhat0 = nan(nx+np,1);
    zhat0(1:nx) = [1e-3*ones(3,1); phat.Base*phat.k_ma/phat.k_tr*ones(4,1); phat.Base];
    zhat0(nx+1:nx+np) = [phat.Base; phat.k_tr; phat.slope; phat.gamma];
    
    % Initialize UKF (log domain)
    jointUkf.State = log(zhat0);
    jointUkf.StateCovariance = blkdiag( diag(ukfInitialStateSDs.^2), diag(ukfInitialParamSDs.^2) );
    
    % 3. Setup GA Optimizer
    ga_options = optimoptions('ga', ...
        'PopulationSize', 30, ...
        'MaxGenerations', 30, ...
        'MaxStallGenerations', 15, ...
        'FunctionTolerance', 1e-3, ...
        'UseParallel', false, ... % IMPORTANT: Disable nested parallelism inside parfor
        'Display', 'off', ...
        'CreationFcn', @createCustomPopulation);
        
    % 4. Closed-loop simulation
    xtrue = x0;
    dose_prev = 0;
    
    hist = struct();
    hist.x = nan(Nk,nx);
    hist.zhat = nan(Nk,nx+np);
    hist.P = cell(Nk,1);
    hist.y = nan(Nk,1);
    hist.u = nan(Nk,1);
    hist.t = Ts_wk*(0:Nk-1);
    hist.patientID = patientIDs(patientIdx);
    hist.trueParams = p;
    
    rng(patientIdx); % Reproducible randomness
    
    for k = 1:Nk
        hist.x(k,:) = xtrue;
        
        % Measurement
        yk = exp( log(xtrue(nx)) + ancMeasNoiseSD*randn );
        hist.y(k) = yk;
        
        % UKF Correction
        try
            correct(jointUkf, log(yk));
        catch ME
            fprintf("UKF correct error (Patient %s, k=%d): %s\n", hist.patientID, k, ME.message);
            break
        end
        
        % Get current estimates
        zhat_current = jointUkf.State;
        P_current = jointUkf.StateCovariance;
        hist.zhat(k,:) = zhat_current;
        hist.P{k} = P_current;
        
        % --- Robust MPC Controller ---
        % if yk < ancSafeLowerLim
        %     uk = 0; % Safety override
        %     disp("ANC below safe limit in iteration " + k);
        % else
            % Define the robust fitness function (Expected cost via UT over parameters)
            % Pass 'phat' (containing fixed parameters) to the fitness function.
            fitnessFcn = @(dose_indices) ga_fitness_robust_params(dose_indices, zhat_current, P_current, phat, Ts_wk, ...
                possible_doses, predictionHorizon, ancWeight, doseWeight, ancTarget, nx, np, ut_params);
            
            try
                % Run GA
                [best_dose_indices, ~] = ga(fitnessFcn, predictionHorizon, ...
                    [], [], [], [], ...
                    ones(1, predictionHorizon), ... % Lower bound
                    n_possible_doses*ones(1, predictionHorizon), ... % Upper bound
                    [], 1:predictionHorizon, ga_options); % Integer constraints
                
                % Implement the first control action
                uk = possible_doses(best_dose_indices(1));
                
            catch ME
                fprintf("GA optimization error (Patient %s, k=%d): %s\n", hist.patientID, k, ME.message);
                uk = dose_prev; % Fallback
            end
        % end
        
        % Apply control and simulate true dynamics
        hist.u(k) = uk;
        try
            xtrue = jost_integrator(xtrue,uk,p,Ts_wk);
        catch ME
             fprintf("Simulation error (Patient %s, k=%d): %s\n", hist.patientID, k, ME.message);
             break;
        end

        % UKF Prediction
        try
            predict(jointUkf, uk);
        catch ME
            fprintf("UKF predict error (Patient %s, k=%d): %s\n", hist.patientID, k, ME.message);
            break
        end
        
        dose_prev = uk;
    end
    
    allResults{patientIdx} = hist;
end

simTime = toc;
fprintf('\nAll simulations completed in %.2f minutes.\n', simTime/60);
fprintf('Average time per patient: %.2f seconds.\n', simTime/numPatients);

%% Save results
saveFilename = sprintf("results/RobustMPC_UT_JostModel_%s.mat", string(datetime("now"),"yyyyMMdd_HHmmss"));
save(saveFilename, 'allResults', 'ancTarget', 'Ts_wk', 'Nk', 'nx', 'np', 'numPatients');
fprintf('Results saved to %s\n', saveFilename);

%% (Optional) Analysis and Plotting

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
% Patients reaching target range (approx. within 30% of target ANC)
targetTolerance = 0.3 * ancTarget;
patientsOnTarget = sum(abs(finalANC - ancTarget) <= targetTolerance);
fprintf('Patients within 30%% of target ANC: %d/%d (%.1f%%)\n', ...
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
        failedSimIdx = [i];
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
    plot(xs,ys, 'LineWidth', 1.5, ...
      'Color', [colors(i,:) 0.4], ...
      'DisplayName', sprintf('Patient %s', string(allResults{i}.patientID)));
end
% xlim([0 max(allResults{1}.t)])
xlabel('Day');
ylabel('6-MP dose (mg/day)');
title('Drug Doses for All Patients (GA Controller)');
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
    plot(xs,ys, 'LineWidth', 1.5, ...
      'Color', [colors(i,:) 0.4], ...
      'DisplayName', sprintf('Patient %s', string(allResults{i}.patientID)));
    
    % Plot measured values:
    plot(allResults{i}.t, allResults{i}.y, ...
      'Color', [colors(i,:) 0.4], 'LineStyle', 'None', 'Marker', 'o', ...
      'DisplayName', sprintf('Patient %s', string(allResults{i}.patientID)));
end
yline(ancTarget, 'b--', 'Target', 'LineWidth', 2);
% xlim([0 max(allResults{1}.t)])
% ylim([0 3]);
xlabel('Day');
ylabel('ANC (10^9 cells/L)');
title('ANC Response for All Patients (GA Controller)');
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
    % ylabel(sprintf('log(x_%d)', j));/
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
paramEstErrors = struct();
finalErrors = struct();
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
    
end

%% Plot relative parameter estimation errors over time
figure('Position', [300 300 1200 800]);
for paramIdx = 1:np
    subplot(2,2,paramIdx);
    hold on;
    
    % Plot error trajectories for all patients
    for i = 1:numPatients
        if isempty(allResults{i}) || size(paramEstErrors.(errorFields{paramIdx}), 1) < i
            continue
        end
        plot(allResults{i}.t, paramEstErrors.(errorFields{paramIdx})(i, :), ...
             'LineWidth', 1, ...
             'DisplayName', sprintf('Patient %s', string(allResults{i}.patientID)));
    end
    
    % Add zero error reference line
    yline(0, 'k--', 'Perfect Estimation', 'LineWidth', 2);
    
    % xlim([0 max(allResults{1}.t)])
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
    
    if isfield(finalErrors, errorFields{paramIdx})
        validErrors = finalErrors.(errorFields{paramIdx})(~isnan(finalErrors.(errorFields{paramIdx})));
    else
        validErrors = [];
    end
    
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

%% ================================================================== %%
%% Helper functions %%
%% ================================================================== %%

% === ROBUST MPC FUNCTIONS ===

% Robust Fitness Function (Expected Cost via UT over Parameter Uncertainty)
function cost = ga_fitness_robust_params(dose_indices, z_mean_initial, P_initial, phat_fixed, Ts_wk, ...
    possible_doses, predictionHorizon, ancWeight, doseWeight, ancTarget, nx, np, ut_params)
    
    % Implements the robust cost function: J = E_theta[ sum (C - y(theta))^2 + R*u^2 ]
    % Expectation is approximated using the Unscented Transform over parameters theta.

    doses = possible_doses(dose_indices);
    
    % 1. Define the uncertainty space
    % State: We condition on the current mean estimate (linear domain).
    x_initial = exp(z_mean_initial(1:nx));
    
    % Parameters: We consider the distribution (log domain).
    p_mean_log = z_mean_initial(nx+1:end);
    P_pp_log = P_initial(nx+1:end, nx+1:end);
    
    % 2. Generate Sigma Points for parameters
    [Xi_p_log, Wmp, ~] = generateSigmaPoints(p_mean_log, P_pp_log, ut_params.Alpha, ut_params.Beta, ut_params.Kappa);
    
    if isempty(Xi_p_log)
        % Sigma point generation failed (e.g., severely ill-conditioned covariance)
        cost = 1e15; % High penalty
        return;
    end
    
    num_sigma_points = size(Xi_p_log, 2);
    expected_cost = 0;

    % 3. Simulate trajectories for each parameter realization (sigma point)
    for j = 1:num_sigma_points
        % Parameter realization (linear domain)
        p_vec_j = exp(Xi_p_log(:, j));
        p_struct_j = makeParamStruct(phat_fixed, p_vec_j);
        
        % Initialize simulation for this realization
        x_current = x_initial;
        cost_j = 0;

        % Simulate over the prediction horizon
        for i = 1:predictionHorizon
            u_step = doses(i);
            
            % Integrate dynamics
            try
                x_next = jost_integrator(x_current, u_step, p_struct_j, Ts_wk);
            catch
                % If simulation fails (e.g., ODE solver issue), assign high cost
                cost_j = cost_j + 1e9; 
                break; % Stop this trajectory
            end

            % Predicted ANC
            anc_predicted = x_next(nx);
            
            % Calculate cost J(theta_j, U)
            cost_j = cost_j + ancWeight * (anc_predicted - ancTarget)^2 ...
                             + doseWeight * (u_step)^2;
            
            x_current = x_next;
        end
        
        % 4. Accumulate expected cost: E[J] = sum_j W_j * J(theta_j, U)
        expected_cost = expected_cost + Wmp(j) * cost_j;
    end
    
    cost = expected_cost;
end

% Generate Sigma Points (Scaled Unscented Transform)
function [Xi, Wm, Wc] = generateSigmaPoints(x_mean, P, alpha, beta, kappa)
    % Generates sigma points and weights for the UT. Includes robust regularization.
    L = length(x_mean);
    lambda = alpha^2 * (L + kappa) - L;
    
    % Ensure P is symmetric
    P = (P + P') / 2;

    % Calculate matrix square root (Cholesky decomposition)
    try
        U = chol((L + lambda) * P, 'lower');
    catch
        % Handle non-positive definite P via regularization (Eigenvalue adjustment)
        [V, D] = eig(P);
        D_diag = diag(D);
        
        % Ensure eigenvalues are non-negative
        min_eig = min(D_diag);
        if min_eig < 1e-9
            % Regularize by shifting eigenvalues slightly positive
            D_diag(D_diag < 1e-9) = 1e-9;
            P_reg = V * diag(D_diag) * V';
        else
            P_reg = P;
        end
        
        try
            U = chol((L + lambda) * P_reg, 'lower');
        catch
            % If regularization fails, return empty (handled by caller)
            Xi = []; Wm = []; Wc = [];
            return;
        end
    end
        
    Xi = zeros(L, 2*L + 1);
    Wm = zeros(1, 2*L + 1);
    Wc = zeros(1, 2*L + 1);
    
    % Sigma point 0 (Mean)
    Xi(:, 1) = x_mean;
    Wm(1) = lambda / (L + lambda);
    Wc(1) = Wm(1) + (1 - alpha^2 + beta);
    
    % Remaining sigma points
    weight_val = 1 / (2 * (L + lambda));
    for i = 1:L
        Xi(:, i+1) = x_mean + U(:, i);
        Xi(:, i+L+1) = x_mean - U(:, i);
        Wm(i+1) = weight_val;
        Wc(i+1) = weight_val;
        Wm(i+L+1) = weight_val;
        Wc(i+L+1) = weight_val;
    end
end


% === MODEL DYNAMICS FUNCTIONS ===

% Integrated Jost 2020 dynamics over the time period Delta
function xnext = jost_integrator(x,u,p,Delta)
    % Integrate PK+PD over Delta (days) with daily bolus dose u.
    if Delta <= 0, xnext = x; return; end
    
    tspan = 0:1:Delta;
    if tspan(end) < Delta, tspan = [tspan Delta]; end
    xnext = x;
    
    % Options for ODE solver. Tighter tolerances than original (1e-1/1e-2) 
    % improve accuracy when exploring diverse parameter spaces.
    options = odeset('RelTol',1e-4, 'AbsTol',1e-5, 'NonNegative', 1:length(x));
    
    deltaC = p.F * u;
    
    for i = 1:length(tspan)-1
        % Apply daily dose
        if mod(tspan(i), 1) < 1e-9 && u > 0
            xnext(1) = xnext(1) + deltaC;
        end
        
        % Integrate over the interval (dose input is zero during integration)
        [~, Y] = ode45(@(t,y) jost_rhs(y, 0, p), [tspan(i) tspan(i+1)], xnext, options);

        if isempty(Y)
             error('ODE solver failed integration.');
        end
        xnext = Y(end,:)';
        xnext = max(xnext, 1e-8); % Ensure positive states
    end
end

% PK/PD dynamics RHS
function dx = jost_rhs(x,dose,p)
    x = max(x,1e-8); % Ensure positive states
    
    % Drug effect and Feedback
    Edrug = p.slope*x(3);
    feedback = (p.Base/x(8))^p.gamma;
    
    dx = zeros(8,1);
    
    % PK (Jost 2020)
    dx(1) = -p.k_a * x(1) + p.F * dose;
    dx(2) = p.k_a * x(1) - p.k_20 * x(2);
    dx(3) = p.FM_3*p.k_me * x(2) - p.CL_6tgn * x(3);
    
    % PD (Friberg/Jost model)
    k_tr = p.k_tr; 
    dx(4) = k_tr*x(4)*(1-Edrug)*feedback - k_tr*x(4);
    dx(5) = k_tr*x(4) - k_tr*x(5);
    dx(6) = k_tr*x(5) - k_tr*x(6);
    dx(7) = k_tr*x(6) - k_tr*x(7);
    dx(8) = k_tr*x(7) - p.k_ma*x(8);
end

% Parameter Structure Helper
function pstruct = makeParamStruct(p_fixed, pvec)
    % Combines fixed parameters (from p_fixed) and variable parameters (from pvec).
    
    pstruct = p_fixed; % Start with fixed parameters
    
    % Overwrite with estimated/sampled parameters (pvec = [Base, k_tr, slope, gamma])
    pstruct.Base = pvec(1);
    pstruct.k_tr = pvec(2);
    pstruct.slope = pvec(3);
    pstruct.gamma = pvec(4);
end

function Population = createCustomPopulation(GenomeLength, FitnessFcn, options)
    % Creates an initial population for the GA.
    % - The first N individuals are constant dose policies, where N is the
    %   number of possible doses.
    % - The remaining individuals are initialized randomly.

    % Get the number of possible doses from the upper bound constraint.
    % Assumes the upper bound is the same for all variables.
    upperBound = options.PopInitRange(2,1);
    lowerBound = options.PopInitRange(1,1);
    
    % Get the total population size.
    populationSize = options.PopulationSize;
    
    % Initialize the population matrix.
    Population = zeros(populationSize, GenomeLength);
    
    % Number of constant-dose individuals to create.
    % This should not exceed the total population size.
    num_const_individuals = min(upperBound, populationSize);
    
    % 1. Create individuals with constant doses.
    % Each individual 'i' will have a dose corresponding to index 'i'.
    for i = 1:num_const_individuals
        Population(i, :) = i;
    end
    
    % 2. Create the rest of the individuals randomly.
    num_random_individuals = populationSize - num_const_individuals;
    if num_random_individuals > 0
        % Generate random integers within the bounds.
        random_part = randi([lowerBound, upperBound], num_random_individuals, GenomeLength);
        
        % Fill the remaining spots in the population.
        Population(num_const_individuals + 1 : end, :) = random_part;
    end
end
