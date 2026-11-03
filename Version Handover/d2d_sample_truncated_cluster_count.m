function n = d2d_sample_truncated_cluster_count(distance_m, model, threshold)
% Sample the number of clusters using a truncated linear Gaussian model.
%
% This function generates the number of clusters (either MPC clusters or
% standalone DMC clusters) as a function of link distance using a linear
% model with Gaussian randomness:
%
%   N(d) = a * d + b + sigma * N(0,1)
%
% The sampled value is continuous and must be converted to a valid integer
% cluster count. A minimum threshold is enforced:
%
%   - threshold = 1 → at least one cluster (MPC clusters)
%   - threshold = 0 → zero or more clusters allowed (standalone DMC)
%
% The procedure:
%   1) Draw a Gaussian realization
%   2) Round to nearest integer
%   3) Reject and re-sample if below threshold
%
% Input:
%   distance_m - Tx–Rx distance (meters)
%   model      - struct containing parameters:
%                a     - slope of cluster count vs distance
%                b     - intercept
%                sigma - standard deviation of Gaussian randomness
%   threshold  - minimum allowed integer value:
%                1 → MPC clusters (n >= 1)
%                0 → standalone DMC clusters (n >= 0)
%
% Output:
%   n          - integer number of clusters satisfying n >= threshold
%
% Example:
%   model = struct('a', 0.2, 'b', 5, 'sigma', 2);
%
%   % MPC clusters (must be >= 1)
%   n_mpc = d2d_sample_truncated_cluster_count(30, model, 1);
%
%   % Standalone DMC clusters (can be 0)
%   n_dmc = d2d_sample_truncated_cluster_count(30, model, 0);

    % ---------------------------------------------------------------------
    % 1) Default threshold (MPC case)
    % ---------------------------------------------------------------------
    if nargin < 3 || isempty(threshold)
        threshold = 1;
    end

    % ---------------------------------------------------------------------
    % 2) Compute mean cluster count based on distance
    % ---------------------------------------------------------------------
    mu = model.a * distance_m + model.b;

    % Standard deviation of variation
    sigma = model.sigma;

    % ---------------------------------------------------------------------
    % 3) Draw samples until threshold condition is satisfied
    % ---------------------------------------------------------------------
    while true
        % Continuous Gaussian realization
        x = mu + sigma * randn;

        % Convert to integer
        n = round(x);

        % Enforce minimum threshold
        if n >= threshold
            return;
        end
    end
end