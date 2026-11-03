function n = d2d_sample_mpc_count(model, threshold)
% Sample the number of MPCs per cluster using a Gamma distribution.
%
% This function generates the number of multipath components (MPCs) within
% a cluster using a Gamma distribution:
%
%   X ~ Gamma(shape = a, scale = b)
%
% The continuous sample is then converted to an integer count. Depending on
% the threshold:
%
%   threshold >= 1:
%       - ensures at least one MPC (standard case)
%       - uses ceil(x) to avoid zero
%
%   threshold = 0:
%       - allows zero MPCs (used for direct cluster extra rays)
%       - uses floor(x) with lower bound at 0
%
% Any sample below the threshold is rejected and re-drawn.
%
% Input:
%   model     - struct with Gamma parameters:
%               a - shape parameter
%               b - scale parameter
%   threshold - minimum allowed integer value:
%               1 -> standard MPC clusters (n >= 1)
%               0 -> direct cluster (n >= 0 allowed)
%
% Output:
%   n         - integer number of MPCs satisfying n >= threshold
%
% Example:
%   model = struct('a', 2, 'b', 3);
%
%   % Standard cluster (at least 1 MPC)
%   n = d2d_sample_mpc_count(model, 1);
%
%   % Direct cluster (allow zero extra MPCs)
%   n = d2d_sample_mpc_count(model, 0);

    % ---------------------------------------------------------------------
    % 1) Default threshold
    % ---------------------------------------------------------------------
    if nargin < 2 || isempty(threshold)
        threshold = 1;
    end

    % ---------------------------------------------------------------------
    % 2) Draw Gamma samples until threshold condition is satisfied
    % ---------------------------------------------------------------------
    while true
        % Continuous Gamma realization
        x = randg(model.a) * model.b;

        % Convert to integer depending on threshold rule
        if threshold <= 0
            % Allow zero (used for direct cluster extra MPCs)
            n = max(0, floor(x));
        else
            % Ensure strictly positive integer
            n = ceil(x);
        end

        % Accept if threshold condition is satisfied
        if n >= threshold
            return;
        end
    end
end