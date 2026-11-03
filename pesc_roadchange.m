
function [out1, out2, out3] = pesc_roadchange(mode, t, z, measurement, p, state)
%PESC_ROADCHANGE  Original PESC plus measurement-triggered road-jump handling.
%
% This NEW file does not replace or modify pesc.m.
%
% Two modes are used by run_pesc_road_change.m:
%
%   [z_new,state_new,detected] = pesc_roadchange('monitor',t,z,meas,p,state)
%      Called once per integration sample to detect an abrupt friction drop
%      based ONLY on changes in measured friction and slip. It resets the
%      PESC estimate and filter states and starts a short freeze interval.
%
%   [zdot,lambda0,xi] = pesc_roadchange('dynamics',t,z,mu,p,state)
%      Uses the exact original pesc.m differential equations except that
%      adaptation and demodulator LPF are held while the freeze is active.
%
% Initialize the additional detector memory in the NEW run script as:
%   state.previous_mu = NaN;
%   state.previous_lambda = NaN;
%   state.freeze_until = -Inf;
%   state.rearm_at = 0;
%   state.detected_times = [];
%
% The measured friction is idealized, just as in the original PESC model.
% The detector is intended for large abrupt friction drops, not gradual roads.

switch mode
    case 'monitor'
        meas = measurement;
        detected = false;
        if isfinite(state.previous_mu)
            change_mu = meas.mu - state.previous_mu;
            change_slip = meas.lambda - state.previous_lambda;
            detected = t >= p.roadctrl.arm_time && t >= state.rearm_at && ...
                change_mu < -p.roadctrl.drop_threshold && ...
                abs(change_slip) <= p.roadctrl.max_lambda_step;
        end
        if detected
            % Fallback is fixed by the NEW controller, not the hidden ice peak.
            z(1) = min(z(1), p.roadctrl.lambda_fallback);
            z(1) = max(0.03, min(0.30, z(1))); % original estimator limits
            z(2) = meas.mu;   % reset HPF state to current measured friction
            z(3) = 0;         % eliminate jump-induced demodulated transient
            state.freeze_until = t + p.roadctrl.freeze_time;
            state.rearm_at = t + p.roadctrl.refractory_time;
            state.detected_times(end+1) = t;
        end
        state.previous_mu = meas.mu;
        state.previous_lambda = meas.lambda;
        out1 = z;
        out2 = state;
        out3 = detected;

    case 'dynamics'
        % Crucial: use UNMODIFIED original PESC with UNMODIFIED parameters.
        [zdot, lambda0, xi] = pesc(t, z, measurement, p);
        if t < state.freeze_until
            zdot(1) = 0;  % freeze slip estimate
            zdot(3) = 0;  % keep LPF gradient reset during freeze
        end
        out1 = zdot;
        out2 = lambda0;
        out3 = xi;

    otherwise
        error('pesc_roadchange:Mode','Unknown mode: %s', mode);
end
end
