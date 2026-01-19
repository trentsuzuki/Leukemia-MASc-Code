function dxdt = JostModel(x,u,popParams,p)
%JOSTMODEL Nonlinear system from Jost, 2020
% x = state
% u = input (6MP dose, mg)
% p = [Base, ktr, slope, gamma]' (patient-specific parameters)

% BSA = 0.82;
% F = 0.22;
% ka = 31.2;
% k20 = 12.72;
% FM3 = 0.019;
% kme = 9.9216;
% CL6tgn = 0.219*BSA^1.16; % From Table 2. Not Equation (2), which is incorrect.
% kma = 2.3765;

Base = p(1);
ktr = p(2);
slope = p(3);
gamma = p(4);

dxdt = zeros(8,1);

dxdt(1) = -popParams.ka*x(1) + popParams.F*u;
dxdt(2) = popParams.ka*x(1) - popParams.k20*x(2);
dxdt(3) = popParams.FM3*popParams.kme*x(2) - popParams.CL6tgn*x(3);
dxdt(4) = ktr*x(4)*(1-slope*x(3))*(Base/x(8))^gamma - ktr*x(4);
dxdt(5) = ktr*(x(4) - x(5));
dxdt(6) = ktr*(x(5) - x(6));
dxdt(7) = ktr*(x(6) - x(7));
dxdt(8) = ktr*x(7) - popParams.kma*x(8);

end

