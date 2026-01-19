%% Generate UKF estimate plots for thesis

clear; clc; close all;

set(0,'defaultTextInterpreter','latex');

% Load results from a previous script run
filename = "ga_ukf_fullJostModel_20250811_0148";
% filename = "ga_ukf_3param_est_20250818_0410";
% filename = "ga_ukf_2param_est_20250819_0943";
% filename = "ga_ukf_fullJostModel_20250827_1021";
load("results/" + filename + ".mat");


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
    % legend('Location', 'eastoutside');
    grid on;
end

%% Calculate and plot relative parameter estimation errors

lineColor = [0.9 0.14 0.14 0.2];
paramNames = {'Base', '$k_{tr}$', 'slope', '$\gamma$'};
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

fig_traj = figure('Visible', 'on');
set(fig_traj, 'Color', 'White', 'Theme', 'Light', 'Position', [100 0 1200 800]);

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
             'color', lineColor, ...
             'DisplayName', sprintf('Patient %s', string(allResults{i}.patientID)));
    end
    
    % Add zero error reference line
    yline(0, 'k--', 'Perfect Estimation', 'LineWidth', 2);
    
    xlim([0 max(allResults{1}.t)])
    xlabel('Day');
    ylabel('Relative Error (\%)');
    title(sprintf('%s', paramNames{paramIdx}), 'fontsize',16);
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
sgtitle("Parameter Estimation Error: $y_r = " + ancTarget + "$ G/L", 'FontSize', 18, 'FontWeight', 'bold');

trajectory_filename = sprintf("figs/" + filename + "_paramEst.png");
exportgraphics(fig_traj, trajectory_filename, 'Resolution', 300);

%% Plot final parameter estimation error distributions

fig_traj = figure('Visible', 'on');
set(fig_traj, 'Color', 'White', 'Theme', 'Light', 'Position', [400 400 1000 600]);

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

%% Helper functions %%

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
    options = odeset('RelTol',1e-4, 'AbsTol',1e-5, 'NonNegative', 1:length(x));
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

function cost = ga_fitness(dose_indices, z_initial, model_fcn, ...
        possible_doses, predictionHorizon, ancWeight, doseWeight, ancTarget, nx)
    doses = possible_doses(dose_indices);
    cost = 0;
    z_current = z_initial;
    
    % Simulate over the prediction horizon
    for i = 1:predictionHorizon
        u_step = doses(i);
        
        % Predict next state using the provided model function
        z_next = model_fcn(z_current, u_step);
        
        % Extract predicted ANC (in log domain)
        anc_predicted = exp(z_next(nx));
        
        % Calculate cost for this step
        cost = cost + ancWeight * (anc_predicted - ancTarget)^2 ...
                    + doseWeight * (u_step)^2;
        
        % Update state for next iteration
        z_current = z_next;
    end
end

function Population = createCustomPopulation(GenomeLength, FitnessFcn, options)
    % createCustomPopulation Creates an initial population for the GA.
    %   - The first N individuals are constant dose policies, where N is the
    %     number of possible doses.
    %   - The remaining individuals are initialized randomly.

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