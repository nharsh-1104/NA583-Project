function [xdot, meas] = abs_plant(x, u, p)
%ABS_PLANT  Single-wheel (quarter-car) braking model, Zhang & Ordonez Eq. (6.1)-(6.3).
%
%   [xdot, meas] = abs_plant(x, u, p)
%
%   States  x = [v; w; s]
%       v : vehicle (wheel-centre) longitudinal speed      [m/s]
%       w : wheel angular speed                            [rad/s]
%       s : distance travelled (for stopping distance)     [m]
%   Input   u : brake torque                               [N m]
%   Params  p : struct from the simulation script (m, R, I, B, g, mu_fun)
%
%   Outputs
%       xdot : state derivative
%       meas : struct of "measured" signals used by the controller / ES
%              .lambda  wheel slip           (v - R w)/v
%              .mu      friction coefficient (what ES maximises)
%              .ax      vehicle acceleration (accelerometer, = dv/dt)
%
%   Note: meas does not depend on u, so the simulation can call this once
%   with u = 0 to get the measurements before computing the control.

v = x(1);
w = x(2);

N      = p.m * p.g;                 % normal load
lambda = (v - p.R * w) / v;         % Eq. (6.3), braking slip in [0, 1]
mu     = p.mu_fun(lambda);          % unknown friction curve

vdot = -N * mu / p.m;                               % Eq. (6.1)
wdot = (-p.B * w + N * p.R * mu - u) / p.I;         % Eq. (6.2)
sdot = v;

xdot = [vdot; wdot; sdot];

meas.lambda = lambda;
meas.mu     = mu;
meas.ax     = vdot;
end
