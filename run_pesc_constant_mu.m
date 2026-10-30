%RUN_PESC_CONSTANT_MU  Perturbation ES for ABS on a constant road surface.
%   Reproduces Zhang & Ordonez (2012), Sec. 6.2, Fig. 6.4.
%   Needs: abs_plant.m, slip_controller.m, pesc.m in the same folder.

clear; clc; close all;

%% ---------------- Parameters ----------------
% Wheel / vehicle (book values; wheel inertia I is not given in the book)
p.m    = 400;          % quarter-car mass             [kg]
p.g    = 9.81;         % gravity                      [m/s^2]
p.R    = 0.3;          % wheel radius                 [m]
p.I    = 1.0;          % wheel inertia (ASSUMED)      [kg m^2]
p.B    = 0.01;         % bearing friction coefficient [N m s]
p.Tmax = 5000;         % actuator torque limit        [N m]

% Constant road: friction curve of Eq. (6.7), peak mu_star at lambda_star
p.lambda_star = 0.05;
p.mu_star     = 0.6;
p.mu_fun = @(lam) 2 * p.mu_star * p.lambda_star * lam ./ (p.lambda_star^2 + lam.^2);

% Inner slip loop
p.c = 20;              % slip-loop bandwidth [rad/s]

% Extremum seeking (book values)
p.es.a   = 0.02;       % dither amplitude
p.es.w   = 10;          % dither frequency     [rad/s]
p.es.wh  = 2;        % high-pass corner     [rad/s]
p.es.wl  = 3;        % low-pass corner      [rad/s] (0 = no LPF)
p.es.k   = 8;        % integrator gain
p.es.phi = 0;          % demodulation phase   [rad]

% Simulation
v0         = 33.33;    % initial speed [m/s] (120 km/h)
lambda_hat0 = 0.1;     % initial optimal-slip guess
v_stop     = 1.0;      % stop simulation below this speed [m/s]
t_end      = 20;       % max simulation time [s]

%% ---------------- Initial conditions ----------------
x0 = [v0; v0/p.R; 0];                 % zero initial slip
z0 = [lambda_hat0; p.mu_fun(lambda_hat0); 0];     % ES states
X0 = [x0; z0];

%% ---------------- Simulate ----------------
opts = odeset('Events', @(t, X) stop_event(t, X, v_stop), ...
              'RelTol', 1e-6, 'AbsTol', 1e-8, 'MaxStep', 1e-3);
[t, X] = ode45(@(t, X) closed_loop(t, X, p), [0 t_end], X0, opts);

%% ---------------- Recompute signals for plotting ----------------
n = numel(t);
[lambda, mu, u, lambda0, lambda_hat, xi] = deal(zeros(n, 1));
for i = 1:n
    [~, sig] = closed_loop(t(i), X(i, :)', p);
    lambda(i)     = sig.lambda;
    mu(i)         = sig.mu;
    u(i)          = sig.u;
    lambda0(i)    = sig.lambda0;
    lambda_hat(i) = X(i, 4);
    xi(i)         = sig.xi;
end
v = X(:, 1);  w = X(:, 2);  s = X(:, 3);

fprintf('Stopping time      : %.2f s (to %.1f m/s)\n', t(end), v_stop);
fprintf('Stopping distance  : %.1f m\n', s(end));
fprintf('Final lambda_hat   : %.3f (true peak %.3f)\n', lambda_hat(end), p.lambda_star);
fprintf('Mean mu (t > 2 s)  : %.3f (peak %.3f)\n', mean(mu(t > 2)), p.mu_star);

%% ---------------- Plots ----------------
figure('Name', 'Perturbation ES - constant mu', 'Position', [100 100 1000 700]);

subplot(3, 2, 1);
plot(t, mu, 'LineWidth', 1.2); hold on;
yline(p.mu_star, '--k', '\mu^*');
xlabel('t [s]'); ylabel('\mu'); title('Friction coefficient'); grid on;

subplot(3, 2, 2);
plot(t, lambda, 'LineWidth', 1.2); hold on;
plot(t, lambda_hat, 'LineWidth', 1.2);
yline(p.lambda_star, '--k', '\lambda^*');
xlabel('t [s]'); ylabel('slip'); legend('\lambda', '\lambda\^_0', 'Location', 'best');
title('Slip and optimal-slip estimate'); grid on;

subplot(3, 2, 3);
plot(t, v, 'LineWidth', 1.2); hold on;
plot(t, p.R * w, 'LineWidth', 1.2);
xlabel('t [s]'); ylabel('[m/s]'); legend('v', 'R\omega', 'Location', 'best');
title('Vehicle and wheel speed'); grid on;

subplot(3, 2, 4);
plot(t, u, 'LineWidth', 1.2);
xlabel('t [s]'); ylabel('T_b [N m]'); title('Brake torque'); grid on;

subplot(3, 2, 5);
plot(t, xi, 'LineWidth', 1.2);
xlabel('t [s]'); ylabel('\xi'); title('Gradient estimate'); grid on;

subplot(3, 2, 6);
lam_grid = linspace(0, 1, 500);
plot(lam_grid, p.mu_fun(lam_grid), 'k', 'LineWidth', 1.2); hold on;
plot(lambda, mu, '.', 'MarkerSize', 4);
xlabel('\lambda'); ylabel('\mu'); title('Operating points on \mu(\lambda)');
legend('\mu(\lambda)', 'trajectory', 'Location', 'best'); grid on;

%% ================= Local functions =================
function [dX, sig] = closed_loop(t, X, p)
    x = X(1:3);     % plant states
    z = X(4:6);     % ES states

    % Measurements (independent of u)
    [~, meas] = abs_plant(x, 0, p);

    % ES: uses measured mu as the performance output
    [zdot, lambda0, xi] = pesc(t, z, meas.mu, p);

    % Inner slip controller tracks lambda0
    u = slip_controller(x, meas, lambda0, p);

    % Plant
    xdot = abs_plant(x, u, p);

    dX = [xdot; zdot];

    sig.lambda  = meas.lambda;
    sig.mu      = meas.mu;
    sig.u       = u;
    sig.lambda0 = lambda0;
    sig.xi      = xi;
end

function [value, isterminal, direction] = stop_event(~, X, v_stop)
    value      = X(1) - v_stop;   % stop when vehicle speed drops to v_stop
    isterminal = 1;
    direction  = -1;
end
