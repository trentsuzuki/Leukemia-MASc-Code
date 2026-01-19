%% Evaluation of model fits from nlmixr2

clear; clc; close all;

if ispc
    userdir = getenv('USERPROFILE') + "/";
else
    userdir = getenv('HOME') + "/";
end

codeFolder = userdir + "UofT/Research/Leukemia/Code/";

% Path to CSV file with patient parameters (Data Sheet 1 from Jost 2020 supplementary materials)
csvFilePath = codeFolder + "Leukemia-PKPD/results/SimplifiedJostModel_20250722_180831/SimplifiedJostModel_20250722_180831.csv";

%% Read patient parameters from CSV
patientData = readtable(csvFilePath);
numPatients = height(patientData);

fprintf('Found %d patients in CSV file\n', numPatients);

%% Read treatment data from CSV

treatmentCsvFilePath = codeFolder + "Jost 2020 Data Analysis/Supplementary Materials/Data Sheet 2.CSV";
treatmentData = readtable(treatmentCsvFilePath);

% For each patient, extract their treatment schedule (days, doses)
patientTreatment = cell(numPatients, 1);
patientNk = zeros(numPatients, 1);
for i = 1:numPatients
    pid = patientData.ID(i);
    mask = treatmentData.ID == pid & treatmentData.EVID ~= 0;
    patientTreatment{i} = treatmentData(mask, :);
    patientNk(i) = height(patientTreatment{i});
end

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

ancMeasNoiseSD = 0;

%% Simulation parameters
nx = 6;
ny = 1;
nu = 1;
np = 4;

ancSafeLowerLim = 0.5; % If ANC is below this level, give zero dose

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
    % Set patient-specific *true* parameters from CSV
    p = p_fixed;
    p.Base = patientBase(patientIdx);
    p.k_tr = patientKtr(patientIdx);
    p.slope = patientSlope(patientIdx);
    p.gamma = patientGamma(patientIdx);
    
    % Initialization
    x0 = [0.1; p.Base*p.k_ma/p.k_tr*ones(4,1); p.Base];
    
    zhat0 = nan(nx+np,1);
    zhat0(1:nx) = [0.1; phat.Base*phat.k_ma/phat.k_tr*ones(4,1); phat.Base];
    zhat0(nx+1:nx+np) = [phat.Base; phat.k_tr; phat.slope; phat.gamma];
    jointUkf.State = log(zhat0);
    jointUkf.StateCovariance = blkdiag( diag(ukfInitialStateSDs.^2), diag(ukfInitialParamSDs.^2) );
    
    % Closed-loop simulation
    xk = x0;
    dose_prev = 0;
    
    hist = struct();
    hist.x = nan(Nk,nx);
    hist.y = nan(Nk,1);
    hist.u = nan(Nk,1);
    hist.t = 0:Nk-1;
    hist.patientID = patientIDs(patientIdx);
    
    for k = 1:Nk
        hist.x(k,:) = xk;
        
        % Measurement (clinic visit)
        hist.y(k) = xk(nx);
        
        % Get dose from treatment history
        uk = patientTreatment{patientIdx}.AMT(k); % Dose for this visit
        hist.u(k) = uk;
        
        % Simulate true dynamics until dose changes
        % T = time until next clinic visit (days)
        xk = friberg_integrator(xk,uk,p,T)
    end
    
    % Store results for this patient
    allResults{patientIdx} = hist;
end

simTime = toc; % End timing
fprintf('\nAll simulations completed in %.2f minutes!\n', simTime/60);
fprintf('Average time per patient: %.2f seconds\n', simTime/numPatients);

%% Save results
save("results/joint_ukf_" + string(datetime("now"),"yyyyMMdd_hhmm"))

%% (Optional) Load results from a previous run

clear;

files = dir("results/joint_ukf_*.mat");
if ~isempty(files)
    [~, idx] = max([files.datenum]);
    latestFile = fullfile("results", files(idx).name);
    fprintf('Loading most recent results: %s\n', latestFile);
    load(latestFile);
else
    fprintf('No results files found in results directory\n');
end

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

function pstruct = makeParamStruct(p,pvec)
    pstruct.ke = p.ke;
    pstruct.F_over_V = p.F_over_V;
    pstruct.Base = pvec(1);
    pstruct.k_tr = pvec(2);
    pstruct.slope = pvec(3);
    pstruct.gamma = pvec(4);
    pstruct.k_ma = p.k_ma;
end