function rays = d2d_generate_mpc_rays(nMPC, centroid, intra, powerDecay, totalPowerLin, d_m, c0, doNormalize)
% Generate specular MPC rays within a single cluster.
%
% This function generates multipath components (MPCs) around a given cluster
% centroid. Each ray is assigned:
%   - an excess delay (relative to the cluster centroid)
%   - angular offsets (AoA and AoD)
%   - a power weight following the intra-cluster decay model
%
% The power model is defined in two regimes (pre/post breakpoint), and the
% resulting linear weights can optionally be normalized.
%
% Input:
%   nMPC           - number of MPCs to generate (integer >= 0)
%   centroid       - struct with fields:
%                    delay_m, aoa_rad, aod_rad
%   intra          - struct with intra-cluster parameters:
%                    tauBP, lambdaDelay, lambdaAngle
%   powerDecay     - struct with fields:
%                    pre, post (each contains zeta, lambdaTau, lambdaPhi)
%   totalPowerLin  - total cluster power assigned to these MPCs (linear scale)
%   d_m            - Tx–Rx distance (meters)
%   c0             - speed of light (m/s)
%   doNormalize    - boolean flag:
%                    true  -> normalize weights to sum to 1
%                    false -> keep raw weights (used in direct cluster)
%
% Output:
%   rays           - struct array of MPC rays, each containing:
%                    delay, angles, power, and phase
%
% Example:
%   rays = d2d_generate_mpc_rays(10, centroid, intra, powerDecay, 1e-6, d_m, c0, true);

    % Default: normalize weights unless explicitly disabled
    if nargin < 8 || isempty(doNormalize)
        doNormalize = true;
    end

    % Initialize empty ray array
    rays = repmat(d2d_make_ray(), 0, 1);

    % If no MPCs requested, return empty
    if nMPC <= 0
        return;
    end

    % ---------------------------------------------------------------------
    % 1) Generate intra-cluster excess delays
    % ---------------------------------------------------------------------
    % deltaTau_m >= 0 (by construction of the segmented model)
    deltaTau_m = d2d_segmented_delay_rnd( ...
        nMPC, intra.tauBP, intra.lambdaDelay(1), intra.lambdaDelay(2), 6);

    % Preallocate raw power weights
    w = zeros(nMPC, 1);

    % ---------------------------------------------------------------------
    % 2) Generate each MPC
    % ---------------------------------------------------------------------
    for m = 1:nMPC

        % Determine which regime (pre/post breakpoint) the ray belongs to
        isPre = deltaTau_m(m) <= intra.tauBP;

        if isPre
            lambdaAng = intra.lambdaAngle(1);
            zeta = powerDecay.pre.zeta;
            ltau = powerDecay.pre.lambdaTau;
            lphi = powerDecay.pre.lambdaPhi;
        else
            lambdaAng = intra.lambdaAngle(2);
            zeta = powerDecay.post.zeta;
            ltau = powerDecay.post.lambdaTau;
            lphi = powerDecay.post.lambdaPhi;
        end

        % -----------------------------------------------------------------
        % Angular offsets (truncated Laplacian, ±60 degrees)
        % -----------------------------------------------------------------
        dAoa = d2d_laplace_rnd(0, lambdaAng, 1, pi/3);
        dAod = d2d_laplace_rnd(0, lambdaAng, 1, pi/3);

        % -----------------------------------------------------------------
        % Power model (in dB domain, then converted to linear)
        % -----------------------------------------------------------------
        p_dB = zeta - ltau * deltaTau_m(m) - lphi * (abs(dAoa) + abs(dAod));
        w(m) = 10^(p_dB / 10);

        % -----------------------------------------------------------------
        % Construct ray struct
        % -----------------------------------------------------------------
        r = d2d_make_ray();
        r.kind = 'mpc';
        r.parent_id = NaN;

        % Excess delay relative to LoS (in meters)
        r.delay_excess_m = centroid.delay_m + deltaTau_m(m);

        % Absolute propagation delay (seconds)
        r.delay_abs_s = (d_m + r.delay_excess_m) / c0;

        % Wrap angles to [0, 2π)
        r.aoa_rad = d2d_wrap_2pi(centroid.aoa_rad + dAoa);
        r.aod_rad = d2d_wrap_2pi(centroid.aod_rad + dAod);

        % Random phase
        r.phase_rad = 2 * pi * rand;

        rays(end+1,1) = r; %#ok<AGROW>
    end

    % ---------------------------------------------------------------------
    % 3) Normalize weights if required
    % ---------------------------------------------------------------------
    if doNormalize
        if sum(w) <= 0
            % Fallback: uniform distribution
            w = ones(nMPC, 1) / nMPC;
        else
            w = w / sum(w);
        end
    end

    % ---------------------------------------------------------------------
    % 4) Assign final power to rays
    % ---------------------------------------------------------------------
    for m = 1:numel(rays)
        rays(m).power_lin = totalPowerLin * w(m);
        rays(m).power_dB = 10 * log10(rays(m).power_lin + eps);
    end
end