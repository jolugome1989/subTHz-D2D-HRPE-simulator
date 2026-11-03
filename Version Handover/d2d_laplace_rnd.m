function x = d2d_laplace_rnd(mu, b, n, boundary)
% Sample truncated Laplace (double-exponential) random variables.
%
% This function generates random samples from a Laplace distribution with:
%   - location parameter mu (center)
%   - scale parameter b (spread)
%
% The standard Laplace distribution has PDF:
%   f(x) = (1 / (2b)) * exp(-|x - mu| / b)
%
% Sampling is performed using inverse transform sampling. The generated
% values are then truncated to lie within a symmetric interval around mu:
%   [mu - boundary, mu + boundary]
%
% Any sample falling outside this interval is rejected and re-drawn. This
% truncation is useful for limiting angular spreads (e.g., ±60°) in channel
% modeling.
%
% Input:
%   mu        - location (mean) of the Laplace distribution
%   b         - scale parameter (controls spread)
%   n         - number of samples to generate (default = 1)
%   boundary  - truncation limit (symmetric around mu), i.e.,
%               output satisfies |x - mu| <= boundary (default = pi)
%
% Output:
%   x         - n-by-1 vector of Laplace-distributed random samples
%               truncated to [mu - boundary, mu + boundary]
%
% Example:
%   % Generate 10 angular offsets within ±60 degrees
%   x = d2d_laplace_rnd(0, deg2rad(10), 10, pi/3);

    % ---------------------------------------------------------------------
    % 1) Handle optional inputs
    % ---------------------------------------------------------------------
    if nargin < 3 || isempty(n)
        n = 1;
    end
    if nargin < 4 || isempty(boundary)
        boundary = pi;   % default full angular range
    end

    % Pre-allocate output
    x = zeros(n, 1);

    % ---------------------------------------------------------------------
    % 2) Generate samples using inverse transform sampling
    % ---------------------------------------------------------------------
    for k = 1:n
        while true
            % Draw uniform random variable in (-0.5, 0.5)
            u = rand - 0.5;

            % Inverse CDF of Laplace distribution
            % This maps u to a Laplace-distributed sample
            sample = mu - b * sign(u) * log(1 - 2 * abs(u));

            % -----------------------------------------------------------------
            % 3) Truncation check
            % -----------------------------------------------------------------
            % Accept only if within the specified boundary
            if abs(sample - mu) <= boundary
                x(k) = sample;
                break;
            end
        end
    end
end