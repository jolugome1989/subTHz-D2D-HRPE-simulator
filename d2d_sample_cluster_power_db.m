function p_dB = d2d_sample_cluster_power_db(ratio, model)
% Sample cluster power in dB using a distance-dependent linear-log model.
%
% This function generates the power of a cluster (in dB) based on a
% log-distance model with Gaussian randomness, as described in the paper.
%
% The model is:
%   P_cl,dB = mu(ratio) + sigma(ratio) * N(0,1)
%
% where:
%   mu(ratio)    = alpha + 10 * beta * log10(ratio)
%   sigma(ratio) = rho   + 10 * eta  * log10(ratio)
%
% The parameter "ratio" typically represents a normalized propagation
% distance, e.g.:
%   ratio = (d + excess_delay) / d
%
% This formulation captures:
%   - mean power decay with distance (controlled by alpha, beta)
%   - increasing/decreasing variability with distance (rho, eta)
%
% Input:
%   ratio - normalized distance (dimensionless, ratio >= 1)
%   model - struct containing parameters:
%           alpha - intercept of mean power (dB)
%           beta  - slope of mean power vs log-distance
%           rho   - baseline standard deviation (dB)
%           eta   - slope of standard deviation vs log-distance
%
% Output:
%   p_dB  - sampled cluster power in dB
%
% Example:
%   ratio = 1.5;
%   model = struct('alpha', -20, 'beta', -1.5, 'rho', 6, 'eta', 0.2);
%   p_dB = d2d_sample_cluster_power_db(ratio, model);

    % ---------------------------------------------------------------------
    % 1) Compute mean cluster power (dB)
    % ---------------------------------------------------------------------
    % Linear function of log-distance
    mu = model.alpha + 10 * model.beta * log10(ratio);

    % ---------------------------------------------------------------------
    % 2) Compute distance-dependent standard deviation
    % ---------------------------------------------------------------------
    sigma = model.rho + 10 * model.eta * log10(ratio);

    % Ensure numerical stability (avoid zero or negative variance)
    sigma = max(sigma, 1e-3);

    % ---------------------------------------------------------------------
    % 3) Generate Gaussian random realization
    % ---------------------------------------------------------------------
    p_dB = mu + sigma * randn;
end