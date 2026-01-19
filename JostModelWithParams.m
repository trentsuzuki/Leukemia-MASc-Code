function dxdt = JostModelWithParams(~,x,u,p)
%JOSTMODEL Nonlinear system from Jost, 2020

dxdt = zeros(8,1);

dxdt(1) = -p.k_a*x(1) + p.F*u;
dxdt(2) = p.k_a*x(1) - p.k_20*x(2);
dxdt(3) = p.FM_3*p.k_me*x(2) - p.CL_6tgn*x(3);
dxdt(4) = p.k_tr*x(4)*(1-p.slope*x(3))*(p.Base/x(8))^p.gamma - p.k_tr*x(4);
dxdt(5) = p.k_tr*(x(4) - x(5));
dxdt(6) = p.k_tr*(x(5) - x(6));
dxdt(7) = p.k_tr*(x(6) - x(7));
dxdt(8) = p.k_tr*x(7) - p.k_ma*x(8);

end