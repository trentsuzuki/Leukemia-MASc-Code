function dxdt = JostModelLastFourStates(~,x,u4,p)
%JOSTMODEL Last four states of the nonlinear system from Jost, 2020
% The input x should be 8x1 for conformity with plotting and analysis scripts.
% The first four states will remain unchanged.

dxdt = zeros(8,1);

dxdt(5) = p.k_tr*(u4 - x(5));
dxdt(6) = p.k_tr*(x(5) - x(6));
dxdt(7) = p.k_tr*(x(6) - x(7));
dxdt(8) = p.k_tr*x(7) - p.k_ma*x(8);

end
