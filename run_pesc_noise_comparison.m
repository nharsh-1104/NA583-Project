%RUN_PESC_NOISE_COMPARISON Four-case ABS/PESC measurement-noise experiment.
% No initialization script needs to run first.
% Required in this folder: abs_plant.m, slip_controller.m, pesc.m,
%                         pesc_robust.m.
% A = original/clean; B = original/noisy;
% C = modified/clean; D = modified/same noise as B.
% Noise enters ONLY the performance measurement supplied to PESC.


clear; clc; close all;
script_folder = fileparts(mfilename('fullpath'));
addpath(script_folder);
required = {'abs_plant.m', 'slip_controller.m', 'pesc.m', 'pesc_robust.m'};
for j = 1:numel(required)
    assert(isfile(fullfile(script_folder, required{j})), ...
        'Missing dependency: %s. Keep all files in the same folder.', required{j});
end

%% 1. Plant and original controller parameters (unchanged)
p.m = 400; p.g = 9.81; p.R = 0.3; p.I = 1.0;
p.B = 0.01; p.Tmax = 5000;
p.lambda_star = 0.05; p.mu_star = 0.60;
p.mu_fun = @(lam) 2*p.mu_star*p.lambda_star*lam ./ ...
    (p.lambda_star^2 + lam.^2);
p.c = 20;
p.es.a = 0.02; p.es.w = 10; p.es.wh = 2;
p.es.wl = 3; p.es.k = 8; p.es.phi = 0;

v0 = 33.33;
lambda_hat0 = 0.10;
v_stop = 1.0;
t_end = 20;
evaluation_start = 2.0;
output_dt = 0.002;

%% 2. Noise and robust-controller settings: edit these for your study
noise_sigma = 0.05;  % SD of Gaussian samples, in friction-coefficient units
noise_Ts = 0.02;     % sample interval [s]; linear interpolation between samples
noise_seed = 7;      % change to repeat the experiment with another realization

p.robust.delta = 0.003;  % dead-zone threshold on FILTERED gradient xi
p.robust.sigma = 0.15;   % leakage rate [1/s]
p.robust.lambda_safe = lambda_hat0; % prior guess; NOT the known true optimum
% These are starting settings, not an optimized or guaranteed robust design.
% A large delta can stop useful adaptation; leakage introduces a bias.

assert(noise_sigma >= 0 && noise_Ts > 0);
assert(p.robust.delta >= 0 && p.robust.sigma >= 0);
assert(p.robust.lambda_safe >= 0.03 && p.robust.lambda_safe <= 0.30);

%% 3. Generate noise ONCE, before calling the ODE solver
rng(noise_seed, 'twister');
noise_time = (0:noise_Ts:(t_end + noise_Ts))';
noise_samples = noise_sigma * randn(size(noise_time));
noise_fun = griddedInterpolant(noise_time, noise_samples, 'linear', 'nearest');
% The same interpolant is used by B and D: n_B(t) == n_D(t).
% The interpolated waveform is colored, finite-bandwidth noise, not ideal
% continuous-time white noise. Its time RMS differs from sample SD.

%% 4. Run all four independent braking simulations
case_names = {'A: Original, clean'; 'B: Original, noisy'; ...
              'C: Modified, clean'; 'D: Modified, noisy'};
modified = [false false true true];
noisy = [false true false true];
x0 = [v0; v0/p.R; 0];
z0 = [lambda_hat0; p.mu_fun(lambda_hat0); 0];
X0 = [x0; z0]; % identical initial states, matching the original script
opts = odeset('Events', @(t,X) stop_event(t,X,v_stop), ...
    'RelTol', 1e-6, 'AbsTol', 1e-8, 'MaxStep', min(1e-3,noise_Ts/5));
results = struct([]);
for j = 1:4
    fprintf('Running %s ...\n', case_names{j});
    rhs = @(t,X) closed_loop(t,X,p,modified(j),noisy(j),noise_fun);
    sol = ode45(rhs, [0 t_end], X0, opts);
    tf = sol.x(end);
    tq = unique([(0:output_dt:tf)'; tf]);
    Xq = deval(sol,tq).';
    signals = zeros(numel(tq),8);
    for i = 1:numel(tq)
        [~,signals(i,:)] = closed_loop(tq(i),Xq(i,:)',p, ...
            modified(j),noisy(j),noise_fun);
    end
    r.name = case_names{j};
    r.t = tq; r.X = Xq;
    r.lambda = signals(:,1); r.mu_true = signals(:,2);
    r.mu_measured = signals(:,3); r.torque = signals(:,4);
    r.lambda0 = signals(:,5); r.xi = signals(:,6);
    r.noise = signals(:,7); r.xi_used = signals(:,8);
    r.lambda_hat = Xq(:,4);
    r.reached_threshold = ~isempty(sol.xe);
    if j == 1
    results = r;       % Initialize using the actual result fields
    else
    results(j) = r;    % Append subsequent cases
    end
    if ~r.reached_threshold
        warning('%s did not reach 1 m/s by t_end.',r.name);
    end
end

%% 5. Fair metrics: use the SAME post-transient time window for all cases
common_end = min(arrayfun(@(r) r.t(end),results));
assert(common_end > evaluation_start, 'Reduce evaluation_start.');
metric_t = unique([(evaluation_start:output_dt:common_end)'; common_end]);
duration = common_end - evaluation_start;
metric_values = zeros(4,8);
reached_threshold = false(4,1);
for j = 1:4
    r = results(j);
    mu = interp1(r.t,r.mu_true,metric_t);
    lh = interp1(r.t,r.lambda_hat,metric_t);
    lam = interp1(r.t,r.lambda,metric_t);
    ref = interp1(r.t,r.lambda0,metric_t);
    tq = interp1(r.t,r.torque,metric_t);
    mean_mu = trapz(metric_t,mu)/duration;
    estimate_rmse = sqrt(trapz(metric_t,(lh-p.lambda_star).^2)/duration);
    tracking_rmse = sqrt(trapz(metric_t,(lam-ref).^2)/duration);
    torque_TV_per_s = sum(abs(diff(tq)))/duration;
    reached_threshold(j) = r.reached_threshold;
    time_to_threshold = NaN; distance_to_threshold = NaN;
    if r.reached_threshold
        time_to_threshold = r.t(end);
        distance_to_threshold = r.X(end,3);
    end
    metric_values(j,:) = [time_to_threshold, distance_to_threshold, ...
        mean_mu, 100*mean_mu/p.mu_star, estimate_rmse, tracking_rmse, ...
        torque_TV_per_s, r.lambda_hat(end)];
end
summary_table = array2table(metric_values,'VariableNames', ...
    {'TimeTo1mps_s','DistanceTo1mps_m','MeanTrueMu','FrictionUse_pct', ...
     'EstimateRMSE','SlipTrackingRMSE','TorqueTVperSecond_Nm_s','FinalSlipEstimate'});
summary_table = addvars(summary_table,case_names,'Before',1,'NewVariableNames','Case');
summary_table.Reached1mps = reached_threshold;
fprintf('\nCommon metric window: %.2f to %.3f s\n',evaluation_start,common_end);
fprintf('Noise sample SD = %.3f; interval = %.3f s; seed = %d\n', ...
    noise_sigma,noise_Ts,noise_seed);
disp(summary_table);
fprintf('B versus D: compare true friction, distance, and estimate RMSE.\n');
fprintf('A versus C: check the clean-performance cost of the modification.\n');
fprintf('One seed is an example, not evidence of general robustness.\n');

%% 6. Plots: identical colors identify cases in every figure
colors = [0 0.35 0.8; 0.85 0.2 0.1; 0.1 0.6 0.3; 0.55 0.15 0.75];
figure('Name','Four-case PESC comparison','Color','w','Position',[60 60 1250 800]);
subplot(3,2,1); plot_cases(results,'mu_true',colors);
yline(p.mu_star,'k--','Peak = 0.60','HandleVisibility','off');
ylabel('\mu_{true}'); title('True friction: performance achieved');
legend(case_names,'Location','best');
subplot(3,2,2); plot_cases(results,'lambda_hat',colors);
yline(p.lambda_star,'k--','Optimum = 0.05','HandleVisibility','off');
ylabel('\hat{\lambda}'); title('Estimated optimal slip');
subplot(3,2,3); plot_cases(results,'lambda',colors);
yline(p.lambda_star,'k--','Optimum = 0.05','HandleVisibility','off');
ylabel('\lambda'); title('Actual wheel slip, including dither response');
subplot(3,2,4); plot_cases(results,'torque',colors);
ylabel('Brake torque [N m]'); title('Brake torque');
subplot(3,2,5); plot_cases(results,'xi',colors);
yline(p.robust.delta,'k--','+dead-zone','HandleVisibility','off');
yline(-p.robust.delta,'k--','-dead-zone','HandleVisibility','off');
ylabel('\xi'); title('Filtered gradient before dead-zone');
subplot(3,2,6); hold on; grid on;
for j = 1:4
    plot(results(j).t,results(j).X(:,1),'Color',colors(j,:),'LineWidth',1.8);
end
xlabel('Time [s]'); ylabel('Vehicle speed [m/s]'); title('Vehicle deceleration');

figure('Name','Noise input and measured friction','Color','w','Position',[100 100 1100 700]);
subplot(3,1,1);
plot(noise_time,noise_samples,'Color',[0.35 0.35 0.35],'LineWidth',1.0);
xlim([0 max(arrayfun(@(r) r.t(end),results))]); grid on;
ylabel('n_\mu'); title('One pre-generated noise waveform shared by B and D');
for k = 1:2
    j = 2*k; r = results(j);
    subplot(3,1,k+1);
    plot(r.t,r.mu_measured,'Color',[0.75 0.75 0.75]); hold on;
    plot(r.t,r.mu_true,'Color',colors(j,:),'LineWidth',1.8);
    yline(p.mu_star,'k--','HandleVisibility','off');
    grid on; ylabel('\mu'); title(r.name);
    legend('Measured \mu supplied to PESC','True \mu','Location','best');
end
xlabel('Time [s]');

figure('Name','Comparison metrics','Color','w','Position',[140 140 1100 600]);
columns = [2 4 5 7];
labels = {'Distance to 1 m/s [m]','Mean friction utilization [%]', ...
          'Optimal-slip estimate RMSE','Torque variation per second [N m/s]'};
for k = 1:4
    subplot(2,2,k);
    b = bar(metric_values(:,columns(k)),'FaceColor','flat'); b.CData = colors;
    set(gca,'XTick',1:4,'XTickLabel',{'A','B','C','D'});
    ylabel(labels{k}); grid on;
    if k == 2, title('Higher is better'); else, title('Lower is better'); end
end

% Results remain in the MATLAB workspace. Optional file export:
% writetable(summary_table, fullfile(script_folder,'noise_comparison_metrics.csv'));
% save(fullfile(script_folder,'noise_comparison_results.mat'), ...
%     'results','summary_table','p','noise_time','noise_samples','noise_seed');

%% Local helpers: MATLAB calls these automatically
function [dX,sig] = closed_loop(t,X,p,modified,noisy,noise_fun)
    x = X(1:3); z = X(4:6);
    [~,meas] = abs_plant(x,0,p);
    n = 0;
    if noisy, n = noise_fun(t); end
    y = meas.mu + n;
    if modified
        [zdot,lambda0,xi] = pesc_robust(t,z,y,p);
        xi_used = sign(xi)*max(abs(xi)-p.robust.delta,0);
    else
        [zdot,lambda0,xi] = pesc(t,z,y,p);
        xi_used = xi;
    end
    % meas.ax stays clean: isolate noise in the outer optimization loop.
    u = slip_controller(x,meas,lambda0,p);
    xdot = abs_plant(x,u,p);
    dX = [xdot;zdot];
    sig = [meas.lambda,meas.mu,y,u,lambda0,xi,n,xi_used];
end

function [value,isterminal,direction] = stop_event(~,X,v_stop)
    value = X(1)-v_stop; isterminal = 1; direction = -1;
end

function plot_cases(results,field,colors)
    hold on; grid on;
    for j = 1:4
        plot(results(j).t,results(j).(field),'Color',colors(j,:),'LineWidth',1.8);
    end
    xlabel('Time [s]');
end
