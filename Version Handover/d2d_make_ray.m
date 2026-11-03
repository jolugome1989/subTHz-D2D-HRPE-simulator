function r = d2d_make_ray()
% Create an empty ray record (struct template for MPC/direct paths).
%
% This function initializes a ray structure with all required fields used
% in the simulator. It serves as a consistent template for storing the
% parameters of specular multipath components (MPCs) and the direct path.
% Each field is initialized with a placeholder value (NaN or empty), which
% will later be populated during the channel generation process.
%
% Input:
%   (none)
%
% Output:
%   r : struct containing the following fields
%       kind            - string indicating ray type ('mpc', 'direct', etc.)
%       cluster_id      - index of the cluster this ray belongs to
%       parent_id       - parent identifier (used for hierarchy if needed)
%       delay_excess_m  - excess delay in meters (relative to LoS path)
%       delay_abs_s     - absolute propagation delay in seconds
%       aoa_rad         - angle of arrival (AoA) in radians
%       aod_rad         - angle of departure (AoD) in radians
%       power_lin       - linear-scale power of the ray
%       power_dB        - power in dB scale
%       phase_rad       - random phase of the ray in radians
%
% Example:
%   r = d2d_make_ray();
%   r.kind = 'mpc';
%   r.power_lin = 1e-6;
%   r.aoa_rad = pi/4;

    r = struct( ...
        'kind', '', ...
        'cluster_id', NaN, ...
        'parent_id', NaN, ...
        'delay_excess_m', NaN, ...
        'delay_abs_s', NaN, ...
        'aoa_rad', NaN, ...
        'aod_rad', NaN, ...
        'power_lin', NaN, ...
        'power_dB', NaN, ...
        'phase_rad', NaN);
end