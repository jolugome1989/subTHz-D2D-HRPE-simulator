function tau = d2d_segmented_delay_rnd(n, bp, lambda1, lambda2, xmax)
% Sample excess delays from a segmented exponential distribution.
%
% This function generates random excess delays according to a two-segment
% (piecewise) exponential model:
%   - Segment 1: [0, bp]        with decay rate lambda1
%   - Segment 2: (bp, xmax]     with decay rate lambda2
%
% The PDF is continuous at the breakpoint (bp) but uses different decay
% rates before and after the breakpoint, reflecting different propagation
% mechanisms within a cluster (e.g., early strong reflections vs. late
% diffuse tails).
%
% The distribution is truncated to [0, xmax]. Any generated sample outside
% this range is rejected and re-drawn.
%
% Input:
%   n        - number of samples to generate
%   bp       - breakpoint (meters) separating the two regimes
%   lambda1  - exponential decay rate for [0, bp]
%   lambda2  - exponential decay rate for (bp, xmax]
%   xmax     - maximum allowed excess delay (meters)
%
% Output:
%   tau      - n-by-1 vector of sampled excess delays (meters), satisfying
%              0 <= tau <= xmax
%
% Example:
%   tau = d2d_segmented_delay_rnd(100, 0.3, 10, 0.5, 6);

    % ---------------------------------------------------------------------
    % 1) Compute normalization constants
    % ---------------------------------------------------------------------
    % Integrals of the two exponential segments:
    %   int1 = ∫_0^bp exp(-lambda1 * x) dx
    %   int2 = ∫_bp^xmax exp(-lambda2 * (x - bp)) dx
    %
    % These are used to normalize the piecewise PDF.
    int1 = (1 - exp(-lambda1 * bp)) / max(lambda1, eps);
    int2 = (1 - exp(-lambda2 * (xmax - bp))) / max(lambda2, eps);

    % Total normalization constant
    C = 1 / (int1 + int2);

    % CDF value at the breakpoint
    % This determines how much probability mass lies in the first segment
    Fbp = C * int1;

    % Preallocate output
    tau = zeros(n, 1);

    % ---------------------------------------------------------------------
    % 2) Generate samples using inverse transform sampling
    % ---------------------------------------------------------------------
    for k = 1:n
        while true
            % Draw a uniform random variable
            u = rand;

            if u <= Fbp
                % ---------------------------------------------------------
                % Segment 1: [0, bp]
                % Inverse CDF of exponential distribution
                % ---------------------------------------------------------
                x = -log(1 - u * lambda1 / C) / lambda1;
            else
                % ---------------------------------------------------------
                % Segment 2: (bp, xmax]
                % Shifted exponential starting at bp
                % ---------------------------------------------------------
                x = bp - log(1 - (u - Fbp) * lambda2 / C) / lambda2;
            end

            % -------------------------------------------------------------
            % 3) Truncation check
            % -------------------------------------------------------------
            % Reject and regenerate if outside valid support
            if x >= 0 && x <= xmax
                tau(k) = x;
                break;
            end
        end
    end
end