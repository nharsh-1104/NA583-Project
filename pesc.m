function [zdot, lambda0, xi] = pesc(t, z, y, p)
%PESC  Perturbation-based extremum seeking controller (Zhang & Ordonez Fig. 6.3).
%
%   [zdot, lambda0, xi] = pesc(t, z, y, p)
%
%   ES states  z = [lambda_hat; y_lp; xi]
%       lambda_hat : current estimate of the optimal slip (integrator state)
%       y_lp       : state of the high-pass filter (a low-pass copy of y)
%       xi         : low-pass filtered demodulated signal ~ (a/2) dmu/dlambda
%
%   Input   y : measured performance output (here mu)
%   Output  lambda0 : slip set-point = lambda_hat + a sin(w t)
%
%   Loop:  y -> HPF s/(s+wh) -> x sin(w t - phi) -> LPF wl/(s+wl) -> k/s -> lambda_hat
%
%   Parameters (in p.es): a, w, wh, wl, k, phi

es = p.es;
lambda_hat = z(1);
y_lp       = z(2);
xi         = z(3);

% 1) High-pass filter: remove the slowly varying mean of y
eta     = y - y_lp;                 % HPF output
y_lp_d  = es.wh * (y - y_lp);       % HPF internal state

% 2) Demodulate with the dither (phi = phase compensation, 0 in the book)
demod = eta * sin(es.w * t - es.phi);

% 3) Low-pass filter (wl = 0 means no filter, as allowed in the book)
if es.wl > 0
    xi_d = es.wl * (demod - xi);
else
    xi   = demod;
    xi_d = 0;
end

% 4) Integrator = gradient ascent on mu
lambda_hat_d = es.k * xi;
if (lambda_hat <= 0.03 && lambda_hat_d < 0) || (lambda_hat >= 0.30 && lambda_hat_d > 0)
    lambda_hat_d = 0;
end

zdot = [lambda_hat_d; y_lp_d; xi_d];

% 5) Add dither to form the slip set-point
lambda0 = lambda_hat + es.a * sin(es.w * t);
end
