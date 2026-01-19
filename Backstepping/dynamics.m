function [dxdt] = dynamics(x,u)
%DYNAMICS Summary of this function goes here
%   Detailed explanation goes here

dxdt = zeros(5,1);

dxdt(1) = -x(1) + u;
dxdt(2) = x(1) - x(2);
dxdt(3) = x(3)*(1-x(2))*nonlinearity(x(5)) - x(3);
dxdt(4) = x(3) - x(4);
dxdt(5) = x(4) - x(5);

end

