function u = slip_controller(x, meas, lambda0, p)
%SLIP_CONTROLLER  Feedback-linearising slip controller, Zhang & Ordonez Eq. (6.5).
%
%   u = slip_controller(x, meas, lambda0, p)
%
%   Makes the closed-loop slip dynamics   d(lambda)/dt = -c (lambda - lambda0)
%   without knowing mu(lambda), because the measured acceleration ax
%   carries the friction information.
%
%   Inputs
%       x       : plant state [v; w; s]
%       meas    : measurements from abs_plant (.lambda, .ax)
%       lambda0 : slip set-point commanded by the ES loop
%       p       : parameters (c, I, R, B, m, Tmax)

v  = x(1);
w  = x(2);
ax = meas.ax;

u = -(p.c * p.I * v / p.R) * (meas.lambda - lambda0) ...
    - p.B * w ...
    - (p.I * w / v) * ax ...
    - p.m * p.R * ax;

% Brake torque cannot be negative (no traction) and is limited by the actuator.
u = min(max(u, 0), p.Tmax);
end
