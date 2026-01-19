% This function defines the system of ordinary differential equations (ODEs)
% for the Friberg model, incorporating the user-specified control law.

function dxdt = friberg_model(~, x, u, params)
    % Unpack state variables
    P = x(1);
    T1 = x(2);
    T2 = x(3);
    T3 = x(4);
    y = x(5); % y is the same as N

    % Unpack parameters
    k_tr = params.k_tr;
    k_ma = params.k_ma;
    Base = params.Base;
    gamma = params.gamma;

    % --- System Dynamics ---
    % The ODEs for the Friberg model under the control law, u
    
    % Progenitor cells dynamics
    dPdt = k_tr * P * ((1-u)*(Base/y)^gamma - 1);
    
    % Transit compartments dynamics
    dT1dt = k_tr * (P - T1);
    dT2dt = k_tr * (T1 - T2);
    dT3dt = k_tr * (T2 - T3);
    
    % Neutrophil (output) dynamics
    dydt = k_tr * T3 - k_ma * y;
    
    % Return the derivatives as a column vector
    dxdt = [dPdt; dT1dt; dT2dt; dT3dt; dydt];
end
