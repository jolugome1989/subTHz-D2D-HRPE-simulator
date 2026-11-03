function cl = d2d_make_dmc_cluster_init()
% Create a template (empty) DMC cluster structure.
%
% This function initializes a DMC cluster struct with all required fields,
% but without assigning any values. It is primarily used for:
%   - preallocation of struct arrays
%   - ensuring consistent field definitions across all DMC clusters
%
% In MATLAB, struct arrays must have identical fields across all elements.
% Therefore, using this template prevents errors such as:
%   "Subscripted assignment between dissimilar structures."
%
% -------------------------------------------------------------------------
% 📤 OUTPUT
% -------------------------------------------------------------------------
% cl : struct with predefined fields (all empty or default values)
%
% Fields:
%   id          → cluster index
%   kind        → 'associated' or 'standalone'
%   associated  → boolean flag
%   parent_id   → ID of parent MPC cluster (if applicable)
%
%   delay_m     → cluster centroid delay (meters)
%   aoa_rad     → angle of arrival (radians)
%   aod_rad     → angle of departure (radians)
%
%   powerLin    → cluster power (linear scale)
%   power_dB    → cluster power (dB scale)
%   sdr_dB      → specular-to-diffuse ratio
%
%   delayDecay  → delay-domain spread parameter
%   lambdaR     → angular spread (Rx)
%   lambdaT     → angular spread (Tx)
%
%   tensor      → DMC tensor [M_f × N_Rx × N_Tx]
%
% -------------------------------------------------------------------------
% 💡 USAGE
% -------------------------------------------------------------------------
% Typical usage for preallocation:
%
%   dmcClusters = repmat(d2d_make_dmc_cluster_init(), 0, 1);
%
% Then populate elements using:
%
%   dmcClusters(end+1) = d2d_make_dmc_cluster(...);
%
% -------------------------------------------------------------------------
% ⚠️ IMPORTANT
% -------------------------------------------------------------------------
% This function ONLY defines structure layout.
% It does NOT assign any physical or simulation values.
%
% All fields must remain consistent with d2d_make_dmc_cluster().
%
% -------------------------------------------------------------------------
% Example:
%   cl = d2d_make_dmc_cluster_init();
% -------------------------------------------------------------------------

    cl = struct( ...
        'id', [], ...
        'kind', '', ...
        'associated', false, ...
        'parent_id', [], ...
        'delay_m', [], ...
        'aoa_rad', [], ...
        'aod_rad', [], ...
        'powerLin', [], ...
        'power_dB', [], ...
        'sdr_dB', [], ...
        'delayDecay', [], ...
        'lambdaR', [], ...
        'lambdaT', [], ...
        'tensor', []);
end