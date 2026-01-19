%% Nonlinear model predictive control using Particle Filter and Genetic Algorithm
% Model from Jost 2020, with NLMPC replaced by a Genetic Algorithm
% State estimation is performed using a Particle Filter (PF) operating in the log-domain.

clear; clc; close all;
userdir = string(java.lang.System.getProperty('user.home')) + "/";
codeFolder = userdir + "UofT/Research/Leukemia/Code/";

% Path to CSV file with patient parameters (Data Sheet 1 from Jost 2020 supplementary materials)
paramFilePath = codeFolder + "Jost 2020 Data Analysis/Supplementary Materials/Data Sheet 1.CSV";
dataFilePath = codeFolder + "Jost 2020 Data Analysis/Supplementary Materials/Data Sheet 2.CSV";

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
p_fixed.k_me = 9.9216; % Metabolic transformation rate constant of 6MP -> 6-TGN
p_fixed.CL_6tgn = nan; % Overwritten in main loop
p_fixed.k_ma = 2.3765; % Death rate of matured neutrophils

% Parameter estimates (same as p_fixed)
phat_fixed.F = p_fixed.F;
phat_fixed.k_a = p_fixed.k_a;
phat_fixed.k_20 = p_fixed.k_20;
phat_fixed.FM_3 = p_fixed.FM_3;
phat_fixed.k_me = p_fixed.k_me;
phat_fixed.CL_6tgn = p_fixed.CL_6tgn;
phat_fixed.k_ma = p_fixed.k_ma;

%% Simulation parameters
Ts_wk = 7; % days
nx = 8;
ny = 1;
nu = 1;
np = 4;
Nk = 30; % weeks to simulate

% Controller parameters
ancTarget = 1.5; % target ANC
ancSafeLowerLim = 0.2; % If ANC is below this level, give zero dose
ancMeasNoiseSD = 0.1; % Standard deviation of measurement noise in log domain
maxAllowableDose = 200; % Maximum 6-MP dose (mg/day)
ancWeight = 1e2; % Weight for ANC
doseWeight = 2e-4; % Weight for 6-MP dose
predictionHorizon = 6; % Number of weeks to predict and control

% Define the quantized dose levels (multiples of 25mg)
possible_doses = 0:25:maxAllowableDose;
n_possible_doses = length(possible_doses);

% Particle Filter parameters (Important: these are for the log-domain state)
nParticles = 5000; % Number of particles for the filter
pfStateProcessNoiseCovs = [1e-6 1e-6 1e-6 5e-4 5e-4 5e-4 5e-4 1e-6];
pfParamProcessNoiseCovs = [1e-4 1e-4 1e-4 1e-4]; % [Base, k_tr, slope, gamma]
pfInitialStateSDs = [1e-2*ones(1,3) 0.125304565*ones(1,4) 0.103687761];
pfInitialParamSDs = 2.*[0.103687761 0.059575249 0.161059312 0.044632609]; % SD of log-params from Jost Data Sheet 1

%% Initialize storage for all patients
allResults = cell(numPatients, 1);

%% Extract patient data for parallel loop (parfor requires simple indexing)
patientIDs = patientData.ID;
patientBase = patientData.Base;
patientKtr = patientData.ktr;
patientSlope = patientData.slope;
patientGamma = patientData.gamma;
% patientBSA is already extracted above

% Leave-one-out averages (calculated in log domain, converted back to linear)
sum_log_Base = sum(log(patientBase));
sum_log_k_tr = sum(log(patientKtr));
sum_log_slope = sum(log(patientSlope));
sum_log_gamma = sum(log(patientGamma));
avg_Base_excl = exp( (sum_log_Base - log(patientBase)) ./ (numPatients - 1) );
avg_ktr_excl = exp( (sum_log_k_tr - log(patientKtr)) ./ (numPatients - 1) );
avg_slope_excl = exp( (sum_log_slope - log(patientSlope)) ./ (numPatients - 1) );
avg_gamma_excl = exp( (sum_log_gamma - log(patientGamma)) ./ (numPatients - 1) );

%% Run simulation for each patient in parallel
fprintf('\nStarting parallel simulation of %d patients...\n', numPatients);
tic; % Start timing

parfor patientIdx = 1:10%%numPatients
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
    
    phat = phat_fixed;
    % Set patient-specific parameter *estimates*
    phat.Base = avg_Base_excl(patientIdx);
    phat.k_tr = avg_ktr_excl(patientIdx);
    phat.slope = avg_slope_excl(patientIdx);
    phat.gamma = avg_gamma_excl(patientIdx);
    % CL_6tgn is considered known
    phat.CL_6tgn = 0.219 * patientBSA(patientIdx)^1.16;
    
    % State transition function for the controller's predictions
    % Input z is the log of the state vector [x; p]
    % This is the deterministic part of the model evolution
    stateFcn = @(z,u) [
        log( jost_integrator(exp(z(1:nx)),u,makeParamStruct(phat,exp(z(nx+1:nx+np))),Ts_wk) );
        z(nx+1:end) % parameters are assumed constant over prediction horizon
    ];
    
    % --- Particle Filter Setup ---
    % Process noise covariance matrix (for log-domain state)
    pfProcessNoiseCov = blkdiag(diag(pfStateProcessNoiseCovs), diag(pfParamProcessNoiseCovs));
    
    % State transition function for the particle filter (operates in log-domain)
    % It propagates each log-particle through the model and adds process noise.
    pfStateTransition = @(prevParticles_log, u) pf_transition_function_log(prevParticles_log, u, stateFcn, pfProcessNoiseCov, nx, np);
    
    % Measurement likelihood function (operates in log-domain)
    % It computes the likelihood of the log-measurement for each predicted log-particle.
    pfMeasurementLikelihood = @(predParticles_log, yk) ...
        normpdf(log(yk), predParticles_log(nx,:), ancMeasNoiseSD);
    
    % Create the particle filter object
    pf = particleFilter(pfStateTransition, pfMeasurementLikelihood);
    pf.StateEstimationMethod = 'maxweight'; %% Set state estimation method
    
    % Initialization

    % True state
    x0 = [1e-3*ones(3,1); p.Base*p.k_ma/p.k_tr*ones(4,1); p.Base]; % 8 elements
    
    % Estimated state
    zhat0 = nan(nx+np,1);
    zhat0(1:nx) = [1e-3*ones(3,1); phat.Base*phat.k_ma/phat.k_tr*ones(4,1); phat.Base];
    zhat0(nx+1:nx+np) = [phat.Base; phat.k_tr; phat.slope; phat.gamma];
    
    % Initialize the filter with particles sampled from the initial log-normal distribution
    log_zhat0 = log(zhat0);
    initialCov = blkdiag(diag(pfInitialStateSDs.^2), diag(pfInitialParamSDs.^2));
    initialize(pf, nParticles, log_zhat0, initialCov); % Initialize with log-domain particles
    
    % Set resampling policy for robustness
    pf.ResamplingPolicy.TriggerMethod = 'interval';
    pf.ResamplingPolicy.SamplingInterval = 1; % Resample every time step
    pf.ResamplingMethod = 'systematic';
    
    % --- GA Options ---
    ga_options = optimoptions('ga', ...
        'PopulationSize', 30, ...
        'MaxGenerations', 30, ...
        'MaxStallGenerations', 15, ...
        'FunctionTolerance', 1e-4, ...
        'UseParallel', false, ...
        'Display', 'off');
        
    % Closed-loop simulation
    xtrue = x0;
    dose_prev = 0;
    
    hist = struct();
    hist.x = nan(Nk,nx);
    hist.zhat = nan(Nk,nx+np);
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
        
        % Particle Filter update with measurement
        try
            correct(pf, yk); % yk is in linear scale, likelihood function handles log transform
        catch ME
            fprintf("Error in PF correct step for patient %s at k=%d: %s\n", hist.patientID, k, ME.message);
            break
        end
        
        % Get the state estimate (mean of log-particles)
        zhat_current_log = getStateEstimate(pf);
        hist.zhat(k,:) = zhat_current_log; % Store log-domain estimate directly
        
        % --- Genetic Algorithm Controller ---
        if yk < ancSafeLowerLim
            uk = 0; % Safety override
        else
            % Define the fitness function for the GA
            % It calculates the cost of a given dose sequence
            fitnessFcn = @(dose_indices) ga_fitness(dose_indices, zhat_current_log, stateFcn, ...
                    possible_doses, predictionHorizon, ancWeight, doseWeight, ancTarget, nx);
            
            try
                % Run the GA to find the best sequence of dose INDICES
                % [best_dose_indices, ~] = ga(fitnessFcn, predictionHorizon, ...
                %     [], [], [], [], ...
                %     ones(1, predictionHorizon), ... % Lower bound (index 1)
                %     n_possible_doses*ones(1, predictionHorizon), ... % Upper bound (index for max dose)
                %     [], 1:predictionHorizon, ... % Integer constraints on all variables
                %     ga_options);
                
                % The control action is the FIRST dose in the optimal sequence
                uk = 100; %% possible_doses(best_dose_indices(1));
                
            catch ME
                fprintf("Error in GA optimization for patient %s at k=%d: %s\n", hist.patientID, k, ME.message);
                uk = dose_prev; % Fallback to previous dose
            end
        end
        
        % Store & apply weekly constant dose
        hist.u(k) = uk;
        
        % Simulate true dynamics over next week
        xtrue = jost_integrator(xtrue,uk,p,Ts_wk);
        
        % Particle Filter time-update for the whole week (Ts_wk) with same uk
        try
            predict(pf, uk);
        catch ME
            fprintf("Error in PF predict step for patient %s at k=%d: %s\n", hist.patientID, k, ME.message);
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
save("results/ga_pf_fullJostModel_" + string(datetime("now"),"yyyyMMdd_hhmm"))

%% (Optional) Load results from a previous run
clear;
files = dir("results/ga_pf_*.mat");
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

%% Helper functions

% Particle Filter State Transition Function (Log Domain)
function new_particles_log = pf_transition_function_log(prev_particles_log, u, stateFcn, Q, nx, np)
    numP = size(prev_particles_log, 2);
    new_particles_log = zeros(size(prev_particles_log));
    
    % Generate all noise samples at once for efficiency
    noise_samples = mvnrnd(zeros(1, nx+np), Q, numP)';
    
    for p_idx = 1:numP
        % Propagate using the deterministic model (in log domain)
        log_z_propagated = stateFcn(prev_particles_log(:, p_idx), u);
        % Add process noise in the log domain
        new_particles_log(:, p_idx) = log_z_propagated + noise_samples(:, p_idx);
    end
end

% Integrated Jost 2020 dynamics over the time period Delta
function xnext = jost_integrator(x,u,p,Delta)
    % Integrate PK+PD over Delta (days) with daily bolus dose u.
    if Delta <= 0
        xnext = x;
        return;
    end
    % Define time points for daily dosing (e.g., 0, 1, 2, ..., Delta)
    tspan = 0:1:Delta;
    if tspan(end) < Delta
        tspan = [tspan Delta]; % Ensure integration covers the full Delta if not an integer
    end
    xnext = x;
    % Options for ODE solver
    options = odeset('RelTol',1e-3, 'AbsTol',1e-5, 'NonNegative', 1:length(x));
    % Pre-calculate the instantaneous concentration increase per dose
    deltaC = p.F * u;
    for i = 1:length(tspan)-1
        % Apply dose at the start of the interval if it aligns with a new day (approx)
        if mod(tspan(i), 1) < 1e-9 && u > 0
            xnext(1) = xnext(1) + deltaC;
        end
        % Integrate over the interval (drug input is zero during integration)
        [~, Y] = ode45(@(t,y) jost_rhs(y, 0, p), [tspan(i) tspan(i+1)], xnext, options);
        if isempty(Y)
            % Handle potential ODE solver failure
            warning('ODE solver failed integration at time %f', tspan(i));
            xnext = nan(size(x));
            return;
        end
        xnext = Y(end,:)';
        xnext = max(xnext, 1e-8); % Ensure positive states (backup)
    end
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

function pstruct = makeParamStruct(p,pvec)
    % Make a struct with all required parameters for the Jost 2020 model
    pstruct.F = p.F;
    pstruct.k_a = p.k_a;
    pstruct.k_20 = p.k_20;
    pstruct.FM_3 = p.FM_3;
    pstruct.k_me = p.k_me;
    pstruct.CL_6tgn = p.CL_6tgn;
    pstruct.Base = pvec(1);
    pstruct.k_tr = pvec(2);
    pstruct.slope = pvec(3);
    pstruct.gamma = pvec(4);
    pstruct.k_ma = p.k_ma;
end

function cost = ga_fitness(dose_indices, z_initial_log, model_fcn, ...
        possible_doses, predictionHorizon, ancWeight, doseWeight, ancTarget, nx)
    doses = possible_doses(dose_indices);
    cost = 0;
    z_current_log = z_initial_log;
    
    % Simulate over the prediction horizon
    for i = 1:predictionHorizon
        u_step = doses(i);
        
        % Predict next state using the provided model function (operates in log domain)
        z_next_log = model_fcn(z_current_log, u_step);
        
        % Extract predicted ANC (convert from log domain)
        anc_predicted = exp(z_next_log(nx));
        
        % Calculate cost for this step
        cost = cost + ancWeight * (anc_predicted - ancTarget)^2 ...
                    + doseWeight * (u_step)^2;
        
        % Update state for next iteration
        z_current_log = z_next_log;
    end
end
