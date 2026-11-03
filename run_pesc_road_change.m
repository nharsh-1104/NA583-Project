
%RUN_PESC_ROAD_CHANGE  Original PESC vs new PESC on sudden dry-to-ice road.
%
% IMPORTANT:
%   This is a NEW simulation file. It does NOT edit the four originals.
%   Original files: abs_plant.m, pesc.m, slip_controller.m,
%                   run_pesc_constant_mu.m
%   New files:      abs_plant_roadchange.m, pesc_roadchange.m,
%                   run_pesc_road_change.m
%
% THREE CASES:
%   1. Original PESC on the original constant surface (validation).
%   2. Original PESC on dry->ice (same original PESC function and gains).
%   3. Road-change PESC on dry->ice (ONLY jump reset/freeze added).
%
% All cases start with the EXACT original parameters, initial states,
% tire-friction curve, and gains. No shared retuning or startup hold.
% Both comparison controllers use the same RK4 numerical integration.
% The constant-road case is also checked against the original ode45 setup.

clear; clc; close all;

%% 1. EXACT original vehicle, inner-loop, and PESC parameters
p.m = 400;
p.g = 9.81;
p.R = 0.3;
p.I = 1.0;
p.B = 0.01;
p.Tmax = 5000;
p.c = 20;
p.lambda_star = 0.05;
p.mu_star = 0.60;
p.mu_fun = @(lam) 2*p.mu_star*p.lambda_star*lam ./ ...
    (p.lambda_star^2 + lam.^2);
p.es.a = 0.02;
p.es.w = 10;
p.es.wh = 2;
p.es.wl = 3;
p.es.k = 8;
p.es.phi = 0;

p.sim.v0 = 33.33;
p.sim.lambda_hat0 = 0.10;
p.sim.v_stop = 1.0;
p.sim.t_end = 30; % longer ONLY to allow low-friction stopping
p.sim.dt = 0.001; % RK4 applied to BOTH vehicle and PESC states

%% 2. New road surface: ORIGINAL road before the boundary, then ice
p.road.s_switch = 45;             % [m], physical distance travelled
p.road.mu_before = p.mu_star;     % EXACT baseline value: 0.60
p.road.lambda_before = p.lambda_star; % EXACT baseline value: 0.05
p.road.mu_after = 0.16;           % illustrative icy-road maximum friction
p.road.lambda_after = 0.03;       % illustrative icy-road optimal slip

%% 3. Settings applying ONLY to the NEW road-change controller
p.roadctrl.lambda_fallback = 0.045; % conservative slip, NOT true optimum
p.roadctrl.drop_threshold = 0.08;   % minimum measured single-step mu drop
p.roadctrl.max_lambda_step = 0.015; % reject large slip jumps
p.roadctrl.arm_time = 0.50;         % detector arming, NOT PESC startup hold
p.roadctrl.freeze_time = 0.30;      % pause new controller adaptation
p.roadctrl.refractory_time = 0.80;  % detector debounce

%% 4. Three tests with an identical initial state and unmodified baseline
p_dry = p;
p_dry.road.s_switch = Inf; % road never changes in baseline validation
original_dry = simulate_case(p_dry, false);
original_change = simulate_case(p, false);
modified_change = simulate_case(p, true);

%% 5. Independent check against ORIGINAL ode45 (dry road only)
% This calls original abs_plant.m, original pesc.m, original slip_controller.m
% and reproduces the initialization/integration method of the original script.
x0 = [p.sim.v0; p.sim.v0/p.R; 0];
z0 = [p.sim.lambda_hat0; p.mu_fun(p.sim.lambda_hat0); 0];
opts = odeset('Events', @(t,X) stop_at_speed(t,X,p.sim.v_stop), ...
    'RelTol',1e-6,'AbsTol',1e-8,'MaxStep',1e-3);
[t_ode, X_ode] = ode45(@(t,X) original_rhs(t,X,p), ...
    [0 p.sim.t_end], [x0;z0], opts);

fprintf('\n===== BASELINE VERIFICATION =====\n');
fprintf('Original ode45 constant road : %.3f s, %.3f m\n', ...
    t_ode(end), X_ode(end,3));
fprintf('New script RK4 constant road : %.3f s, %.3f m\n', ...
    original_dry.t(end), original_dry.s(end));
fprintf('Difference (RK4 - ode45)    : %.4f m\n', ...
    original_dry.s(end)-X_ode(end,3));

fprintf('\n===== SUDDEN ROAD CHANGE COMPARISON =====\n');
fprintf('Dry  peak: mu*=%.2f, lambda*=%.3f (your ORIGINAL road)\n', ...
    p.road.mu_before, p.road.lambda_before);
fprintf('Ice  peak: mu*=%.2f, lambda*=%.3f\n', ...
    p.road.mu_after, p.road.lambda_after);
fprintf('Road boundary: %.1f m, sample step: %.4f s\n', ...
    p.road.s_switch,p.sim.dt);
print_case('Original PESC: dry only',original_dry);
print_case('Original PESC: dry -> ice',original_change);
print_case('New PESC: dry -> ice',modified_change);
if ~isempty(modified_change.detected_times)
    fprintf('New controller detected at: %.4f s\n', ...
        modified_change.detected_times(1));
else
    fprintf('WARNING: the new controller did not detect the road transition.\n');
end

%% 6. Plots: original dry baseline (gray), original road-change (red),
%%          new road-change controller (blue)
figure('Name','PESC dry-to-ice comparison with original baseline', ...
    'Color','w','Position',[80 70 1250 850]);

subplot(3,2,1);
plot(original_dry.t,original_dry.mu,'Color',[0.45 0.45 0.45]); hold on;
plot(original_change.t,original_change.mu,'r','LineWidth',1.1);
plot(modified_change.t,modified_change.mu,'b','LineWidth',1.1);
plot(original_change.t,original_change.mu_peak_true,'k--');
mark_change(original_change.switch_time);
ylabel('\mu'); xlabel('Time (s)'); title('Friction coefficient'); grid on;
legend('Original: constant road','Original: road change', ...
    'New: road change','Available peak: road change','Location','best');

subplot(3,2,2);
plot(original_dry.t,original_dry.lambda,'Color',[0.45 0.45 0.45]); hold on;
plot(original_change.t,original_change.lambda,'r');
plot(modified_change.t,modified_change.lambda,'b');
plot(original_change.t,original_change.lambda_peak_true,'k--');
mark_change(original_change.switch_time);
ylabel('\lambda'); xlabel('Time (s)'); title('Actual slip'); grid on;
legend('Original: constant road','Original: road change', ...
    'New: road change','True road optimum','Location','best');

subplot(3,2,3);
plot(original_dry.t,original_dry.lambda_hat,'Color',[0.45 0.45 0.45]); hold on;
plot(original_change.t,original_change.lambda_hat,'r');
plot(modified_change.t,modified_change.lambda_hat,'b');
plot(original_change.t,original_change.lambda_peak_true,'k--');
mark_change(original_change.switch_time);
if ~isempty(modified_change.detected_times)
    xline(modified_change.detected_times(1),':b','Detected');
end
ylabel('Estimated optimal slip'); xlabel('Time (s)');
title('PESC optimum estimate'); grid on;
legend('Original: constant road','Original: road change', ...
    'New: road change','True road optimum','Location','best');

subplot(3,2,4);
plot(original_dry.t,original_dry.v,'Color',[0.45 0.45 0.45]); hold on;
plot(original_change.t,original_change.v,'r');
plot(modified_change.t,modified_change.v,'b');
mark_change(original_change.switch_time);
ylabel('Vehicle speed (m/s)'); xlabel('Time (s)');
title('Vehicle speed'); grid on;
legend('Original: constant road','Original: road change', ...
    'New: road change','Location','best');

subplot(3,2,5);
plot(original_dry.t,original_dry.u,'Color',[0.45 0.45 0.45]); hold on;
plot(original_change.t,original_change.u,'r');
plot(modified_change.t,modified_change.u,'b');
mark_change(original_change.switch_time);
ylabel('Brake torque (N m)'); xlabel('Time (s)');
title('Applied brake torque'); grid on;
legend('Original: constant road','Original: road change', ...
    'New: road change','Location','best');

subplot(3,2,6);
plot(original_dry.t,original_dry.s,'Color',[0.45 0.45 0.45]); hold on;
plot(original_change.t,original_change.s,'r');
plot(modified_change.t,modified_change.s,'b');
yline(p.road.s_switch,'k--','Road boundary');
mark_change(original_change.switch_time);
ylabel('Distance (m)'); xlabel('Time (s)');
title('Stopping distance'); grid on;
legend('Original: constant road','Original: road change', ...
    'New: road change','Location','best');

saveas(gcf,'pesc_road_change_comparison.png');
save('pesc_road_change_results.mat','p', 'original_dry', ...
    'original_change','modified_change','t_ode','X_ode');

%% ================= NEW helper functions only =================
function out = simulate_case(p, use_roadchange_controller)
    h = p.sim.dt;
    num_samples = floor(p.sim.t_end/h)+1;
    % CRITICAL: EXACT ORIGINAL initial PESC filter state, not zero!
    X = [p.sim.v0; p.sim.v0/p.R; 0; p.sim.lambda_hat0; ...
        p.mu_fun(p.sim.lambda_hat0); 0];

    state.previous_mu = NaN;
    state.previous_lambda = NaN;
    state.freeze_until = -Inf;
    state.rearm_at = 0;
    state.detected_times = [];

    fields = {'t','v','w','s','lambda','mu','lambda0','lambda_hat', ...
        'xi','u','mu_peak_true','lambda_peak_true'};
    for j = 1:numel(fields)
        out.(fields{j}) = zeros(num_samples,1);
    end
    n = num_samples;
    for i = 1:num_samples
        t = (i-1)*h;
        [~, meas] = abs_plant_roadchange(X(1:3),0,p);
        if use_roadchange_controller
            [z_new,state] = pesc_roadchange('monitor', ...
                t,X(4:6),meas,p,state);
            X(4:6) = z_new;
        end
        [~, signal] = closed_loop(t,X,p,use_roadchange_controller,state);
        out.t(i) = t;
        out.v(i) = X(1);
        out.w(i) = X(2);
        out.s(i) = X(3);
        out.lambda(i) = meas.lambda;
        out.mu(i) = meas.mu;
        out.lambda0(i) = signal.lambda0;
        out.lambda_hat(i) = X(4);
        out.xi(i) = signal.xi;
        out.u(i) = signal.u;
        out.mu_peak_true(i) = meas.mu_peak_true;
        out.lambda_peak_true(i) = meas.lambda_peak_true;

        if X(1) <= p.sim.v_stop
            n = i;
            break;
        end
        % Same RK4 integration of ALL six states in both controller cases.
        a = closed_loop(t,X,p,use_roadchange_controller,state);
        b = closed_loop(t+h/2,X+h*a/2,p,use_roadchange_controller,state);
        c = closed_loop(t+h/2,X+h*b/2,p,use_roadchange_controller,state);
        d = closed_loop(t+h,X+h*c,p,use_roadchange_controller,state);
        X = X + (h/6)*(a+2*b+2*c+d);
    end
    for j = 1:numel(fields)
        vals = out.(fields{j});
        out.(fields{j}) = vals(1:n);
    end
    changed = find(out.s >= p.road.s_switch,1,'first');
    if isempty(changed)
        out.switch_time = NaN;
    else
        out.switch_time = out.t(changed);
    end
    out.detected_times = state.detected_times;
    out.stopped = out.v(end) <= p.sim.v_stop;
end

function [dX, signal] = closed_loop(t,X,p,enhanced,state)
    x = X(1:3);
    z = X(4:6);
    [~,meas] = abs_plant_roadchange(x,0,p);
    if enhanced
        [zdot,lambda0,xi] = pesc_roadchange('dynamics', ...
            t,z,meas.mu,p,state);
    else
        [zdot,lambda0,xi] = pesc(t,z,meas.mu,p); % ORIGINAL, unchanged
    end
    u = slip_controller(x,meas,lambda0,p); % ORIGINAL, unchanged
    xdot = abs_plant_roadchange(x,u,p);
    dX = [xdot;zdot];
    signal.lambda0 = lambda0;
    signal.xi = xi;
    signal.u = u;
end

function dX = original_rhs(t,X,p)
    x = X(1:3);
    z = X(4:6);
    [~,meas] = abs_plant(x,0,p); % ORIGINAL plant
    [zdot,lambda0] = pesc(t,z,meas.mu,p); % ORIGINAL PESC
    u = slip_controller(x,meas,lambda0,p); % ORIGINAL slip loop
    xdot = abs_plant(x,u,p);
    dX = [xdot;zdot];
end

function [value,isterminal,direction] = stop_at_speed(~,X,v_stop)
    value = X(1)-v_stop;
    isterminal = 1;
    direction = -1;
end

function mark_change(time)
    if isfinite(time)
        xline(time,'k:','Road change');
    end
end

function print_case(name,d)
    fprintf('%s\n',name);
    if d.stopped
        fprintf('  Time to 1 m/s     : %.3f s\n',d.t(end));
        fprintf('  Distance to 1 m/s : %.3f m\n',d.s(end));
    else
        fprintf('  Did not reach 1 m/s before simulation ends.\n');
    end
    fprintf('  Final slip estimate: %.4f\n',d.lambda_hat(end));
    if isfinite(d.switch_time)
        fprintf('  Road changed at   : %.3f s\n',d.switch_time);
    end
    fprintf('  Change detections : %d\n',numel(d.detected_times));
end
