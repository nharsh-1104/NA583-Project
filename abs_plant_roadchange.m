function [xdot, meas] = abs_plant_roadchange(x, u, p)
%ABS_PLANT_ROADCHANGE  Quarter-car longitudinal ABS plant with a road jump.
%   Same dynamics as the original abs_plant.m, except that the tire-road
%   friction curve is selected by TRAVELLED DISTANCE, x(3).
%
%   x = [vehicle speed (m/s); wheel speed (rad/s); distance (m)]
%   u = applied brake torque (N*m)
%
%   The road parameters below are PLANT truth; the adaptive controller
%   receives only the measured .mu and .lambda fields.

v = x(1);
w = x(2);
s = x(3);

if s < p.road.s_switch
    mu_peak = p.road.mu_before;
    lambda_peak = p.road.lambda_before;
else
    mu_peak = p.road.mu_after;
    lambda_peak = p.road.lambda_after;
end

v_safe = max(v, 1e-3); % avoid division by zero at very low speed
lambda = (v - p.R*w) / v_safe;

% Same analytic friction curve used in the original constant-road model:
mu = 2*mu_peak*lambda_peak*lambda / (lambda_peak^2 + lambda^2);

normal_load = p.m*p.g;
vdot = -normal_load*mu/p.m;
wdot = (-p.B*w + normal_load*p.R*mu - u)/p.I;
sdot = v;

xdot = [vdot; wdot; sdot];
meas.lambda = lambda;
meas.mu = mu;
meas.ax = vdot;
% Ground-truth entries used ONLY for evaluating the simulation, NOT by PESC.
meas.lambda_peak_true = lambda_peak;
meas.mu_peak_true = mu_peak;
end
