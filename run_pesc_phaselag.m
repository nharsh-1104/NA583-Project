%RUN_PESC_CASE1_PHASELAG  Case 1: Phase Lag Failure Mode & Demodulator Fix
%   Reproduces Zhang & Ordonez (2012) PESC with Actuator Phase Lag.
%   Needs: abs_plant_actuator.m, slip_controller.m, pesc.m

clear; clc; close all;

%% ---------------- Parameters ----------------
p.m    = 400;          % quarter-car mass             [kg]
p.g    = 9.81;         % gravity                      [m/s^2]
p.R    = 0.3;          % wheel radius                 [m]
p.I    = 1.0;          % wheel inertia                [kg m^2]
p.B    = 0.01;         % bearing friction coefficient [N m s]
p.Tmax = 5000;         % actuator torque limit        [N m]

% Actuator lag & inner slip loop dynamics
p.tau_a = 0.08;        % Actuator time constant [s] (causes phase lag)
p.c     = 10;          % Inner slip loop bandwidth [rad/s]

% Friction curve
p.lambda_star = 0.15;  % true peak slip
p.mu_star     = 0.8;   % true peak friction
p.mu_fun = @(lam) 2 * p.mu_star * p.lambda_star * lam ./ (p.lambda_star^2 + lam.^2);

% PESC parameters
p.es.a   = 0.02;       % dither amplitude
p.es.w   = 15;         % dither frequency [rad/s]
p.es.wh  = 2;          % high-pass corner [rad/s]
p.es.wl  = 3;          % low-pass corner  [rad/s]
p.es.k   = 5;          % integrator gain

% Calculate total phase lag phi at dither frequency w
phi_slip = atan(p.es.w / p.c);
phi_act  = atan(p.es.w * p.tau_a);
phi_hpf  = atan(p.es.wh / p.es.w); % phase lead from HPF
phi_calc = phi_slip + phi_act - phi_hpf;

fprintf('Dither frequency w = %.1f rad/s\n', p.es.w);
fprintf('Calculated Phase Lag = %.2f rad (%.1f deg)\n', phi_calc, rad2deg(phi_calc));

% Simulation setup
v0          = 33.33;   % 120 km/h
lambda_hat0 = 0.05;    % initial guess below peak
v_stop      = 1.0;    
t_end       = 10;      

%% ---------------- Initial Conditions ----------------
% State mapping: X = [v; w; s; u; lambda_hat; y_lp; xi]
x0 = [v0; v0/p.R; 0; 0];                      % plant states [v; w; s; u]
z0 = [lambda_hat0; p.mu_fun(lambda_hat0); 0];  % ES states [lambda_hat; y_lp; xi]
X0 = [x0; z0];

opts = odeset('Events', @(t, X) stop_event(t, X, v_stop), 'RelTol', 1e-5, 'MaxStep', 1e-3);

%% ---------------- Run 1: Uncompensated (phi = 0) ----------------
p_un = p;
p_un.es.phi = 0; % no phase compensation
[t_un, X_un] = ode45(@(t, X) closed_loop(t, X, p_un), [0 t_end], X0, opts);

%% ---------------- Run 2: Compensated (phi = phi_calc) ----------------
p_co = p;
p_co.es.phi = phi_calc; % phase compensation applied
[t_co, X_co] = ode45(@(t, X) closed_loop(t, X, p_co), [0 t_end], X0, opts);

%% ---------------- Recompute Output Signals ----------------
% Recompute signals for Uncompensated run
n_un = numel(t_un);
[lam_un, mu_un, u_un, lam0_un, xi_un] = deal(zeros(n_un, 1));
for i = 1:n_un
    [~, sig] = closed_loop(t_un(i), X_un(i, :)', p_un);
    lam_un(i)  = sig.lambda;
    mu_un(i)   = sig.mu;
    u_un(i)    = sig.u;
    lam0_un(i) = sig.lambda0;
    xi_un(i)   = sig.xi;
end

% Recompute signals for Compensated run
n_co = numel(t_co);
[lam_co, mu_co, u_co, lam0_co, xi_co] = deal(zeros(n_co, 1));
for i = 1:n_co
    [~, sig] = closed_loop(t_co(i), X_co(i, :)', p_co);
    lam_co(i)  = sig.lambda;
    mu_co(i)   = sig.mu;
    u_co(i)    = sig.u;
    lam0_co(i) = sig.lambda0;
    xi_co(i)   = sig.xi;
end

%% ---------------- All 6 Subplots Grid ----------------
figure('Name', 'Case 1: Phase Lag Comparison (6 Subplots)', 'Position', [100 80 1100 800]);

% 1. Friction Coefficient
subplot(3, 2, 1);
plot(t_un, mu_un, 'r--', 'LineWidth', 1.2); hold on;
plot(t_co, mu_co, 'b', 'LineWidth', 1.2);
yline(p.mu_star, '--k', '\mu^*');
xlabel('t [s]'); ylabel('\mu'); title('Friction coefficient'); grid on;
legend('Uncomp (\phi=0^o)', sprintf('Comp (\\phi=%.0f^o)', rad2deg(phi_calc)), 'Location', 'best');

% 2. Slip and Optimal Slip Estimate
subplot(3, 2, 2);
plot(t_un, lam_un, 'r:', 'LineWidth', 0.8); hold on;
plot(t_un, X_un(:, 5), 'r--', 'LineWidth', 1.2); % X(:,5) is lambda_hat
plot(t_co, lam_co, 'c:', 'LineWidth', 0.8);
plot(t_co, X_co(:, 5), 'b', 'LineWidth', 1.2);
yline(p.lambda_star, '--k', '\lambda^*');
xlabel('t [s]'); ylabel('slip'); title('Slip and optimal-slip estimate'); grid on;
legend('\lambda (uncomp)', '\lambda_0 (uncomp)', '\lambda (comp)', '\lambda_0 (comp)', 'Location', 'best');

% 3. Vehicle and Wheel Speed
subplot(3, 2, 3);
plot(t_un, X_un(:, 1), 'r--', 'LineWidth', 1.2); hold on;
plot(t_un, p.R * X_un(:, 2), 'm:', 'LineWidth', 1.0);
plot(t_co, X_co(:, 1), 'b', 'LineWidth', 1.2);
plot(t_co, p.R * X_co(:, 2), 'c:', 'LineWidth', 1.0);
xlabel('t [s]'); ylabel('[m/s]'); title('Vehicle and wheel speed'); grid on;
legend('v (uncomp)', 'R\omega (uncomp)', 'v (comp)', 'R\omega (comp)', 'Location', 'best');

% 4. Brake Torque
subplot(3, 2, 4);
plot(t_un, u_un, 'r--', 'LineWidth', 1.2); hold on;
plot(t_co, u_co, 'b', 'LineWidth', 1.2);
xlabel('t [s]'); ylabel('T_b [N m]'); title('Brake torque'); grid on;
legend('Uncomp', 'Comp', 'Location', 'best');

% 5. Gradient Estimate
subplot(3, 2, 5);
plot(t_un, xi_un, 'r--', 'LineWidth', 1.2); hold on;
plot(t_co, xi_co, 'b', 'LineWidth', 1.2);
yline(0, 'k:');
xlabel('t [s]'); ylabel('\xi'); title('Gradient estimate'); grid on;
legend('Uncomp', 'Comp', 'Location', 'best');

% 6. Operating Points on mu(lambda)
subplot(3, 2, 6);
lam_grid = linspace(0, 1, 500);
plot(lam_grid, p.mu_fun(lam_grid), 'k', 'LineWidth', 1.2); hold on;
plot(lam_un, mu_un, 'r.', 'MarkerSize', 3);
plot(lam_co, mu_co, 'b.', 'MarkerSize', 3);
xlabel('\lambda'); ylabel('\mu'); title('Operating points on \mu(\lambda)'); grid on;
legend('\mu(\lambda)', 'Uncomp traj', 'Comp traj', 'Location', 'best');

%% ================= Local functions =================
function [dX, sig] = closed_loop(t, X, p)
    x = X(1:4);     % plant states [v; w; s; u]
    z = X(5:7);     % ES states    [lambda_hat; y_lp; xi]

    % Measurements (plant independent of u_cmd)
    [~, meas] = abs_plant_actuator_lag(x, 0, p);

    % Outer PESC loop
    [zdot, lambda0, xi] = pesc(t, z, meas.mu, p);

    % Inner slip controller calculates u_cmd
    u_cmd = slip_controller(x, meas, lambda0, p);

    % Plant step
    xdot = abs_plant_actuator_lag(x, u_cmd, p);

    dX = [xdot; zdot];

    sig.lambda  = meas.lambda;
    sig.mu      = meas.mu;
    sig.u       = meas.u;
    sig.lambda0 = lambda0;
    sig.xi      = xi;
end

function [val, term, dir] = stop_event(~, X, v_stop)
    val  = X(1) - v_stop; % stop when vehicle speed drops to v_stop
    term = 1;
    dir  = -1;
end