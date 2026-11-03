function sim = subthz_d2dhr_simulator(cfg)
%SUBTHZ_D2DHR_SIMULATOR
% Stochastic HRPE-based sub-THz outdoor D2D channel simulator.
%
% -------------------------------------------------------------------------
% 📖 GENERATION SEQUENCE (CRITICAL)
% -------------------------------------------------------------------------
% 1) Generate MPC rays (including direct cluster if LoS)
% 2) Build MPC DDIR (frequency → delay transform)
% 3) Determine threshold from MPC-only channel
% 4) Generate DMC clusters (threshold-aware)
% 5) Combine MPC and DMC contributions
%
% NOTE:
% The ordering is important:
%   - Threshold must come from MPC (not DMC)
%   - DMC scaling depends on this threshold
%   - Avoids circular dependency between MPC and DMC

%% ------------------------------------------------------------------------
% 0) Initialization
% ------------------------------------------------------------------------
if nargin < 1 || isempty(cfg)
    cfg = struct();
    % define the configuration of the sounder for the emulated channel
    % if the cfg is not empty (has valid input), then properly assign
    % otherwise, take the default values
    % cfg = d2d_apply_defaults(cfg);

    % Fix random seed → ensures reproducibility of stochastic simulation
    % rng(cfg.seed);
    % rng(73);
end

% define the configuration of the sounder for the emulated channel
% if the cfg is not empty (has valid input), then properly assign
% otherwise, take the default values
cfg = d2d_apply_defaults(cfg);

% Fix random seed → ensures reproducibility of stochastic simulation
rng(cfg.seed);

% ParaSounder = d2d_build_parasounder(cfg);
% Warn user if simulation is outside measurement-supported range
if cfg.distance_m < 1 || cfg.distance_m > 100
    warning('Distance outside recommended range [1,100] m.');
end

%% ------------------------------------------------------------------------
% 1) Load parameters
% ------------------------------------------------------------------------
% All statistical parameters (Gamma, Lognormal, etc.) come from measurement
p = d2d_params();

scenarioKey = lower(cfg.scenario);  % 'los' or 'nlos'
if ~ismember(scenarioKey,{'los','nlos'})
    error("cfg.scenario must be 'LOS' or 'NLOS'.")
end

% Tx-Rx separation distance in meters
d_m = cfg.distance_m;

% Speed of light
c0 = 299792458;

% Wavelength
lambda_m = c0 / cfg.fc_Hz;

% Free-space path loss (reference power level)
% Used to anchor all generated cluster powers
fspl_dB = 20 * log10(4*pi*d_m/lambda_m);
pRefLin = 10^(-fspl_dB/10);

%% ------------------------------------------------------------------------
% 2) Grid definition
% ------------------------------------------------------------------------
% Delay grid (in meters, later mapped to seconds)
delayGrid_m = 0:cfg.delayStep_m:cfg.maxDelay_m;

% Angular grid (azimuth only)
phiGridRad = deg2rad(-cfg.dmcAngleMaxDeg:cfg.dmcAngleStepDeg:cfg.dmcAngleMaxDeg);

% Pack grid structure
grid.delay_m = delayGrid_m(:);
grid.phi_R   = phiGridRad(:).';
grid.phi_T   = phiGridRad(:).';
grid.delta_phi = deg2rad(cfg.dmcAngleStepDeg);

% Dimensions
M_f  = numel(grid.delay_m);   % delay bins
N_Rx = cfg.N_Rx;              % Rx antennas
N_Tx = cfg.N_Tx;              % Tx antennas

%% ------------------------------------------------------------------------
% 3) Storage initialization
% ------------------------------------------------------------------------
% Dynamic arrays for rays and clusters
mpcRays = repmat(d2d_make_ray(),0,1);
dmcClusters = repmat(d2d_make_dmc_cluster_init(),0,1);

clusterId = 0;

% Channel tensors
ddir_mpc = complex(zeros(M_f,N_Rx,N_Tx));  % specular
ddir_dmc = complex(zeros(M_f,N_Rx,N_Tx));  % diffuse

%% ------------------------------------------------------------------------
% 4) DIRECT CLUSTER (LoS ONLY)
% ------------------------------------------------------------------------
% This block explicitly models the LoS path + its stochastic perturbations
if strcmpi(scenarioKey,'los')

    clusterId = clusterId + 1;

    % Direct path is centered at zero excess delay and boresight angles
    directCentroid.delay_m = 0;
    directCentroid.aoa_rad = 0;
    directCentroid.aod_rad = 0;

    % Rician fading captures LoS amplitude fluctuations
    directAmp = d2d_sample_rician_amp( ...
        p.direct.rician.nu, ...
        p.direct.rician.sigma);

    % Absolute LoS power
    P_direct_abs = pRefLin * directAmp^2;

    % Additional MPCs around LoS (small perturbations)
    nExtraMPC = d2d_sample_mpc_count(p.mpcCount.direct, 0);

    % Generate them WITHOUT normalization → raw power shares
    extraRays = d2d_generate_mpc_rays( ...
        nExtraMPC, ...
        directCentroid, ...
        p.intra.direct, ...
        p.direct.powerDecay, ...
        1.0, ...
        d_m, ...
        c0, ...
        false);

    % Total stochastic share
    P_stochastic_share = sum([extraRays.power_lin]);

    % Reject unrealistic cases (stochastic dominating LoS)
    while (P_stochastic_share >= 1)
        extraRays = d2d_generate_mpc_rays( ...
            nExtraMPC, directCentroid, p.intra.direct, ...
            p.direct.powerDecay, 1.0, d_m, c0, false);
        P_stochastic_share = sum([extraRays.power_lin]);
    end

    % Residual ensures total MPC power is preserved
    P_direct_norm = max(eps, 1 - P_stochastic_share);

    % Total MPC power of direct cluster
    P_MPC_abs = P_direct_abs / P_direct_norm;

    % Scale stochastic MPCs accordingly
    for k = 1:numel(extraRays)
        extraRays(k).power_lin = extraRays(k).power_lin * P_MPC_abs;
        extraRays(k).power_dB = 10 * log10(extraRays(k).power_lin + eps);
        extraRays(k).cluster_id = clusterId;
        mpcRays(end+1,1) = extraRays(k);
    end

    % Add main LoS ray
    directRay = d2d_make_ray();
    directRay.kind = 'direct';
    directRay.cluster_id = clusterId;
    directRay.delay_excess_m = 0;
    directRay.delay_abs_s = d_m / c0;
    directRay.aoa_rad = 0;
    directRay.aod_rad = 0;
    directRay.power_lin = P_direct_abs;
    directRay.power_dB = 10 * log10(P_direct_abs + eps);
    directRay.phase_rad = 2*pi*rand;

    mpcRays(end+1,1) = directRay;
end

%% ------------------------------------------------------------------------
% 5) MPC CLUSTERS (GENERAL)
% ------------------------------------------------------------------------
% These represent non-LoS specular clusters
mpcCountModel = p.mpcClusterCount.(scenarioKey);

if strcmpi(scenarioKey,'los')
    % LoS already has one cluster → allow zero extra clusters
    nMpcClusters = d2d_sample_truncated_cluster_count(d_m, mpcCountModel, 0);
else
    nMpcClusters = d2d_sample_truncated_cluster_count(d_m, mpcCountModel, 1);
end

for c = 1:nMpcClusters

    clusterId = clusterId + 1;

    % Cluster centroid delay and angles
    centroid.delay_m = d2d_sample_gamma(p.mpcArrival.(scenarioKey).a,...
                                        p.mpcArrival.(scenarioKey).b);

    centroid.aoa_rad = d2d_sample_cluster_angle(p.clusterAngles.(scenarioKey));
    centroid.aod_rad = d2d_sample_cluster_angle(p.clusterAngles.(scenarioKey));

    % Number of rays inside cluster
    nMPC = d2d_sample_mpc_count(p.mpcCount.(scenarioKey),1);

    % Distance-dependent cluster power
    ratio = max((d_m + centroid.delay_m)/d_m,1);
    pc_dB = d2d_sample_cluster_power_db(ratio,p.clusterPower.(scenarioKey));
    pcLin = pRefLin * 10^(pc_dB/10);

    % Generate rays (normalized inside function)
    specRays = d2d_generate_mpc_rays( ...
        nMPC, centroid, p.intra.(scenarioKey), ...
        p.mpcPowerDecay.(scenarioKey), pcLin, d_m, c0, true);

    for k = 1:numel(specRays)
        specRays(k).cluster_id = clusterId;
        mpcRays(end+1,1) = specRays(k);
    end
end

%% ------------------------------------------------------------------------
% 6) Build MPC DDIR tensor
% ------------------------------------------------------------------------
% Convert discrete rays into continuous delay-domain response

delayStep_s = cfg.delayStep_m / c0;

% Frequency resolution consistent with delay sampling
deltaF_Hz = 1/(M_f * delayStep_s);

% Centered frequency grid → avoids time-domain shift
freqGrid_Hz = ((0:M_f-1) - floor(M_f/2)).' * deltaF_Hz;

Hf_mpc = complex(zeros(M_f,N_Rx,N_Tx));

for r = 1:numel(mpcRays)

    % Array response (EADF interpolation)
    rayField = d2d_path_field_from_ray(mpcRays(r),cfg);

    % Convert delay → phase rotation
    tau_s = mpcRays(r).delay_excess_m / c0;

    phaseVec = exp(-1j*2*pi*freqGrid_Hz*tau_s);

    % Outer product across frequency and spatial dimensions
    Hf_mpc = Hf_mpc + reshape(phaseVec,[M_f,1,1]) ...
                        .* reshape(rayField,[1,N_Rx,N_Tx]);
end

% IFFT → delay domain (produces sinc-like spreading)
ddir_mpc = ifft(ifftshift(Hf_mpc,1),[],1);

%% ------------------------------------------------------------------------
% 7) Threshold from MPC
% ------------------------------------------------------------------------
% Maximum observed power defines dynamic range reference
Pmax = max(abs(ddir_mpc(:)).^2);

% Linear threshold
thresholdLin = Pmax * 10^(-cfg.dynamicRange_dB/10);

%% ------------------------------------------------------------------------
% 8) DMC GENERATION
% ------------------------------------------------------------------------
standaloneCountModel = p.standaloneDMCCount.(scenarioKey);
nStandaloneDMC = d2d_sample_truncated_cluster_count(d_m, standaloneCountModel, 0);

for s = 1:nStandaloneDMC

    centroid.delay_m = d2d_sample_gamma(p.standaloneDMC.arrival.(scenarioKey).a,...
                                        p.standaloneDMC.arrival.(scenarioKey).b);

    centroid.aoa_rad = d2d_sample_cluster_angle(p.standaloneDMC.angles.(scenarioKey));
    centroid.aod_rad = d2d_sample_cluster_angle(p.standaloneDMC.angles.(scenarioKey));

    ratio = max((d_m + centroid.delay_m)/d_m,1);
    pdmc_dB = d2d_sample_cluster_power_db(ratio,p.standaloneDMC.power.(scenarioKey));
    pdmcLin = pRefLin * 10^(pdmc_dB/10);

    % Generate correlated diffuse tensor
    [hDMC, info] = d2d_generate_dmc_cluster_tensor( ...
        centroid, p.standaloneDMC.spread.(scenarioKey), ...
        pdmcLin, grid, cfg, thresholdLin);

    ddir_dmc = ddir_dmc + hDMC;

    clusterId = clusterId + 1;

    dmcClusters(end+1,1) = d2d_make_dmc_cluster( ...
        clusterId,'standalone',false,centroid,...
        pdmcLin,NaN,info,hDMC);
end

%% ------------------------------------------------------------------------
% 9) Combine
% ------------------------------------------------------------------------
ddir = ddir_mpc + ddir_dmc;

%% ------------------------------------------------------------------------
% 10) Output
% ------------------------------------------------------------------------
sim.mpcRays = mpcRays;
sim.dmcClusters = dmcClusters;
sim.ddir_mpc = ddir_mpc;
sim.ddir_dmc = ddir_dmc;
sim.ddir = ddir;
sim.grid = grid;
sim.cfg = cfg;

sim.meta.distance_m = d_m;
sim.meta.scenario = upper(cfg.scenario);
sim.meta.dynamicRange_dB = cfg.dynamicRange_dB;
sim.meta.referencePowerLinear = pRefLin;

end