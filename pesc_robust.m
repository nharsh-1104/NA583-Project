function [zdot, lambda0, xi] = pesc_robust(t, z, y, p)
%PESC_ROBUST Original PESC filters/dither with dead-zone and leakage.
% No knowledge of p.lambda_star or p.mu_star is used by this controller.
% Set p.robust.delta = 0 and p.robust.sigma = 0 to recover original PESC.

es = p.es;
lambda_hat = z(1);
y_lp = z(2);
xi = z(3);

eta = y - y_lp;
y_lp_d = es.wh * (y - y_lp);
demod = eta * sin(es.w * t - es.phi);
if es.wl > 0
    xi_d = es.wl * (demod - xi);
else
    xi = demod;
    xi_d = 0;
end

% Continuous dead-zone: remove small updates, not the dither itself.
xi_used = sign(xi) * max(abs(xi) - p.robust.delta, 0);
% Leakage toward an independently selected reference (not the true peak).
lambda_hat_d = es.k * xi_used ...
    - p.robust.sigma * (lambda_hat - p.robust.lambda_safe);

% Same outward-update guard as the original controller.
if (lambda_hat <= 0.03 && lambda_hat_d < 0) || ...
        (lambda_hat >= 0.30 && lambda_hat_d > 0)
    lambda_hat_d = 0;
end

zdot = [lambda_hat_d; y_lp_d; xi_d];
lambda0 = lambda_hat + es.a * sin(es.w * t);
end
