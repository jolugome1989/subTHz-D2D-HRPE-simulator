function cfg = d2d_apply_defaults(cfg)
% Merge user-supplied configuration with default parameters.
%
% This function fills in missing fields in the user-provided configuration
% struct `cfg` using predefined default values, without overwriting fields
% that are already specified and non-empty.
%
% In addition, this function:
%   - supports backward-compatible aliases (N_R, N_T)
%   - prepares internal parameters required by the EADF evaluator
%   - validates key parameters (e.g., dynamic range, EADF dimensions)
%
% Input:
%   cfg - user-defined configuration struct (possibly incomplete)
%
% Output:
%   cfg - completed and validated configuration struct
%
% Example:
%   cfg = struct('distance_m', 50, 'dynamicRange_dB', 40);
%   cfg = d2d_apply_defaults(cfg);
%   % dynamicRange_dB will be clipped to 35 dB with a warning

    % ---------------------------------------------------------------------
    % 1) Define default configuration
    % ---------------------------------------------------------------------
    def = struct();

    % Basic simulation parameters
    def.distance_m = 35;           % Tx–Rx distance (meters)
    def.scenario = 'LOS';          % propagation scenario
    def.fc_Hz = 145.5e9;          % carrier frequency (Hz)
    def.seed = 1;                 % RNG seed

    % Array configuration
    def.N_Rx = 36;                % number of receive antennas
    def.N_Tx = 36;                % number of transmit antennas

    % Delay grid
    def.maxDelay_m = 300;         % maximum delay window (meters)
    def.delayStep_m = 0.3;        % delay resolution (meters)

    % Angular grid for DMC
    def.dmcAngleMaxDeg = 180;     % angular range (degrees)
    def.dmcAngleStepDeg = 10;     % angular resolution (degrees)

    % Dynamic range (important constraint)
    def.dynamicRange_dB = 35;     % maximum allowed dynamic range

    % Overlap rejection thresholds
    def.separationTau_m = 1.0;    % delay separation
    def.separationAngle_deg = 15; % angular separation

    % EADF inputs
    def.EADF_RX = [];
    def.EADF_TX = [];

    % EADF evaluator settings
    def.sampling_mode = 'uniform';
    def.ArrayFrequencyDependence = 'FrequencyIndependent';
    def.estimationmode = 'azimuth';

    % Backward-compatible aliases
    def.N_R = [];
    def.N_T = [];

    % ---------------------------------------------------------------------
    % 2) Fill missing fields (do NOT overwrite user-defined values)
    % ---------------------------------------------------------------------
    fields = fieldnames(def);
    for k = 1:numel(fields)
        f = fields{k};
        if ~isfield(cfg, f) || isempty(cfg.(f))
            cfg.(f) = def.(f);
        end
    end

    % ---------------------------------------------------------------------
    % 3) Handle backward-compatible aliases
    % ---------------------------------------------------------------------
    if isfield(cfg, 'N_R') && ~isempty(cfg.N_R) && ...
       (~isfield(cfg, 'N_Rx') || isempty(cfg.N_Rx))
        cfg.N_Rx = cfg.N_R;
    end

    if isfield(cfg, 'N_T') && ~isempty(cfg.N_T) && ...
       (~isfield(cfg, 'N_Tx') || isempty(cfg.N_Tx))
        cfg.N_Tx = cfg.N_T;
    end

    % Ensure dimensions are defined
    if ~isfield(cfg, 'N_Rx') || isempty(cfg.N_Rx)
        cfg.N_Rx = 36;
    end
    if ~isfield(cfg, 'N_Tx') || isempty(cfg.N_Tx)
        cfg.N_Tx = 36;
    end

    % ---------------------------------------------------------------------
    % 4) Internal aliases required by EADF evaluator
    % ---------------------------------------------------------------------
    cfg.N_R = cfg.N_Rx;
    cfg.N_T = cfg.N_Tx;

    % Number of delay samples
    cfg.M_f = numel(0:cfg.delayStep_m:cfg.maxDelay_m);

    % Disable normalization control inside simulator
    cfg.NormFactorPhi_RX = 1;
    cfg.NormFactorTheta_RX = 1;
    cfg.NormFactorPhi_TX = 1;
    cfg.NormFactorTheta_TX = 1;

    % ---------------------------------------------------------------------
    % 5) Validate dynamic range
    % ---------------------------------------------------------------------
    % The model assumes dynamic range should not exceed 35 dB.
    % If the user provides a larger value, clip it and issue a warning.

    if cfg.dynamicRange_dB > 35
        warning('d2d_apply_defaults:DynamicRangeExceeded', ...
            'dynamicRange_dB (%.2f dB) exceeds 35 dB. Clipping to 35 dB.', ...
            cfg.dynamicRange_dB);
        cfg.dynamicRange_dB = 35;
    end

    % Ensure dynamic range is non-negative
    if cfg.dynamicRange_dB < 0
        warning('d2d_apply_defaults:DynamicRangeInvalid', ...
            'dynamicRange_dB (%.2f dB) is negative. Resetting to 0 dB.', ...
            cfg.dynamicRange_dB);
        cfg.dynamicRange_dB = 0;
    end

    % ---------------------------------------------------------------------
    % 6) Validate EADF matrices
    % ---------------------------------------------------------------------
    d2d_validate_eadf_dimensions(cfg);
end