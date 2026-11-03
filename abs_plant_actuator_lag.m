function [xdot, meas] = abs_plant_actuator_lag(x, u_cmd, p)
%ABS_PLANT_ACTUATOR Plant model with first-order hydraulic brake actuator dynamics.
%   States x = [v; w; s; u] (u = actual brake torque [N m])

v = x(1);
w = x(2);
u = x(4); % actual brake torque state

N      = p.m * p.g;                 
lambda = (v - p.R * w) / v;         
mu     = p.mu_fun(lambda);          

vdot = -N * mu / p.m;                               
wdot = (-p.B * w + N * p.R * mu - u) / p.I;         
sdot = v;
udot = (u_cmd - u) / p.tau_a;        % first-order actuator dynamic

xdot = [vdot; wdot; sdot; udot];

meas.lambda = lambda;
meas.mu     = mu;
meas.ax     = vdot;
meas.u      = u;
end