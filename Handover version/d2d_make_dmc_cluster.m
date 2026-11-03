function cl = d2d_make_dmc_cluster(clusterId, kind, associated, centroid, powerLin, sdr_dB, info, tensor)
% Create a DMC cluster structure (Diffuse Multipath Component).
%
% This function constructs a structured representation of a DMC cluster,
% which models the diffuse (non-specular) scattering component of the
% channel. Each DMC cluster is characterized by:
%   - a centroid in delay and angular domain
%   - total power
%   - spatial spreading properties
%   - a full tensor representation in delay × Rx × Tx domain
%
% The function is designed to produce a consistent struct layout, which is
% critical for building arrays of clusters in MATLAB.
%
% -------------------------------------------------------------------------
% 📥 INPUT
% -------------------------------------------------------------------------
% clusterId   : unique identifier for this cluster
%
% kind        : string describing cluster type
%               'associated' → linked to an MPC cluster
%               'standalone' → independent DMC cluster
%
% associated  : boolean flag
%               true  → associated DMC (shares centroid with MPC cluster)
%               false → standalone DMC
%
% centroid    : struct containing cluster center parameters:
%               delay_m  → excess delay (meters)
%               aoa_rad  → angle of arrival (radians)
%               aod_rad  → angle of departure (radians)
%
% powerLin    : total cluster power (linear scale)
%
% sdr_dB      : specular-to-diffuse ratio in dB
%               (NaN for standalone DMC clusters)
%
% info        : struct containing additional model parameters:
%               delayDecay → delay-domain decay constant
%               lambdaR    → angular spread parameter (Rx)
%               lambdaT    → angular spread parameter (Tx)
%
% tensor      : complex-valued DMC tensor:
%               size → [M_f × N_Rx × N_Tx]
%               represents diffuse channel contribution
%
% -------------------------------------------------------------------------
% 📤 OUTPUT
% -------------------------------------------------------------------------
% cl          : struct containing DMC cluster description:
%
%   id         → cluster index
%   kind       → cluster type ('associated' / 'standalone')
%   associated → logical flag
%   parent_id  → parent MPC cluster ID (0 if standalone)
%
%   delay_m    → centroid delay
%   aoa_rad    → centroid AoA
%   aod_rad    → centroid AoD
%
%   powerLin   → cluster power (linear)
%   power_dB   → cluster power (dB)
%   sdr_dB     → SDR value (if applicable)
%
%   delayDecay → delay spread parameter
%   lambdaR    → angular spread (Rx)
%   lambdaT    → angular spread (Tx)
%
%   tensor     → full DMC contribution tensor
%
% -------------------------------------------------------------------------
% 💡 NOTES
% -------------------------------------------------------------------------
% - Associated DMC clusters share their centroid with MPC clusters
% - Standalone DMC clusters are independently generated
% - power_dB uses log10 with eps to avoid numerical issues
%
% - This function assumes all inputs are valid and consistent
%   (validation should be done earlier in the pipeline)
%
% -------------------------------------------------------------------------
% Example:
%   cl = d2d_make_dmc_cluster(1, 'standalone', false, centroid, 1e-6, NaN, info, tensor);
% -------------------------------------------------------------------------

    % Initialize struct
    cl = struct();

    % ---------------------------------------------------------------------
    % Basic identifiers
    % ---------------------------------------------------------------------
    cl.id = clusterId;
    cl.kind = kind;
    cl.associated = associated;

    % Parent cluster ID (0 means no parent / standalone)
    cl.parent_id = 0;

    % ---------------------------------------------------------------------
    % Centroid parameters
    % ---------------------------------------------------------------------
    cl.delay_m = centroid.delay_m;
    cl.aoa_rad = centroid.aoa_rad;
    cl.aod_rad = centroid.aod_rad;

    % ---------------------------------------------------------------------
    % Power representation
    % ---------------------------------------------------------------------
    cl.powerLin = powerLin;

    % Convert to dB (numerically safe)
    cl.power_dB = 10 * log10(powerLin + eps);

    % Specular-to-diffuse ratio
    cl.sdr_dB = sdr_dB;

    % ---------------------------------------------------------------------
    % Spread parameters
    % ---------------------------------------------------------------------
    cl.delayDecay = info.delayDecay;
    cl.lambdaR = info.lambdaR;
    cl.lambdaT = info.lambdaT;

    % ---------------------------------------------------------------------
    % Tensor representation (core DMC contribution)
    % ---------------------------------------------------------------------
    cl.tensor = tensor;

end