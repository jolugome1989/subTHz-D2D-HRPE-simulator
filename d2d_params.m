function p = d2d_params()
% Return the numerical parameter set used by the simulator.
%
% This function collects the fitted model parameters reported in the paper
% draft. The parameters are grouped by modeling purpose:
%   - number of MPC clusters / standalone DMC clusters
%   - number of MPCs per cluster
%   - cluster arrival-time models
%   - centroid angle models
%   - intra-cluster excess-delay and angular-offset models
%   - direct-path Rician fading
%   - SDR models for cluster power splitting
%   - intra-cluster power decay models
%   - cluster power models
%   - associated and standalone DMC spread models
%
% Output:
%   p : struct containing all simulation parameters used by the channel
%       generator.

    % ---------------------------------------------------------------------
    % 1) Number of clusters versus link distance
    % ---------------------------------------------------------------------
    % The paper models the number of MPC clusters and standalone DMC
    % clusters with distance-dependent truncated linear Gaussian laws.
    % Each entry stores:
    %   a     - slope
    %   b     - intercept
    %   sigma - standard deviation of the Gaussian perturbation
    %
    % The simulator later samples an integer cluster count and rejects values
    % below 1.
    p.mpcClusterCount.los  = struct('a', 0.17,  'b', 10.22, 'sigma', 7.40);
    p.mpcClusterCount.nlos = struct('a', -0.45, 'b', 80.98, 'sigma', 33.63);

    p.standaloneDMCCount.los  = struct('a', 0.21,  'b', 39.56, 'sigma', 17.10);
    p.standaloneDMCCount.nlos = struct('a', -0.82, 'b', 77.83, 'sigma', 31.77);

    % ---------------------------------------------------------------------
    % 2) Number of MPCs per cluster
    % ---------------------------------------------------------------------
    % The number of MPCs inside each cluster is modeled by a Gamma law.
    % These shape/scale parameters are used for:
    %   - direct cluster extra MPCs
    %   - LoS MPC clusters
    %   - NLoS MPC clusters
    %
    % The sampled result is converted to an integer count in the simulator.
    p.mpcCount.direct = struct('a', 3.32, 'b', 3.57);
    p.mpcCount.los    = struct('a', 0.96, 'b', 6.69);
    p.mpcCount.nlos   = struct('a', 0.78, 'b', 24.39);

    % ---------------------------------------------------------------------
    % 3) Cluster arrival-time distributions
    % ---------------------------------------------------------------------
    % Cluster centroids in the delay domain are modeled using Gamma
    % distributions. These are the excess-delay centroids of the MPC clusters
    % and standalone DMC clusters.
    p.mpcArrival.los  = struct('a', 1.47, 'b', 29.85);
    p.mpcArrival.nlos = struct('a', 1.78, 'b', 40.07);

    p.standaloneDMC.arrival.los  = struct('a', 2.78, 'b', 31.32);
    p.standaloneDMC.arrival.nlos = struct('a', 4.58, 'b', 27.29);

    % ---------------------------------------------------------------------
    % 4) Cluster centroid angle distributions
    % ---------------------------------------------------------------------
    % LoS cluster centroids follow a bimodal Laplacian model, representing
    % two dominant angular lobes.
    % NLoS cluster centroids are modeled with a uniform angular distribution.
    %
    % Fields:
    %   mode   - angular model type
    %   w      - mixture weights for the two Laplacian lobes
    %   mu     - lobe centers
    %   lambda - Laplace decay parameters
    p.clusterAngles.los = struct( ...
        'mode', 'bimodal_laplace', ...
        'w', [0.46, 0.54], ...
        'mu', [0, pi], ...
        'lambda', [0.20, 1.14]);

    p.clusterAngles.nlos = struct( ...
        'mode', 'uniform', ...
        'w', [], 'mu', [], 'lambda', []);

    % Standalone DMC cluster centroids use the same angular modeling
    % philosophy, but with their own fitted parameters.
    p.standaloneDMC.angles.los = struct( ...
        'mode', 'bimodal_laplace', ...
        'w', [0.56, 0.44], ...
        'mu', [0, pi], ...
        'lambda', [0.29, 0.56]);

    p.standaloneDMC.angles.nlos = struct( ...
        'mode', 'uniform', ...
        'w', [], 'mu', [], 'lambda', []);

    % ---------------------------------------------------------------------
    % 5) Intra-cluster MPC excess delay and angular offset models
    % ---------------------------------------------------------------------
    % These parameters define the segmented excess-delay law and the angular
    % spread of MPCs within each cluster.
    %
    % Fields:
    %   tauBP       - breakpoint between the pre- and post-breakpoint regimes
    %   lambdaDelay - decay parameters for the two delay segments
    %   lambdaAngle - angular spread parameters for the two regimes
    %
    % The direct cluster uses its own fitted values, and LoS/NLoS clusters use
    % their own separate parameter sets.
    p.intra.direct = struct( ...
        'tauBP', 0.30, ...
        'lambdaDelay', [11.63, 0.37], ...
        'lambdaAngle', [0.14, 0.19]);

    p.intra.los = struct( ...
        'tauBP', 0.37, ...
        'lambdaDelay', [9.36, 0.61], ...
        'lambdaAngle', [0.10, 0.18]);

    p.intra.nlos = struct( ...
        'tauBP', 0.29, ...
        'lambdaDelay', [6.85, 0.33], ...
        'lambdaAngle', [0.15, 0.24]);

    % ---------------------------------------------------------------------
    % 6) Direct-path Rician fading and SDR
    % ---------------------------------------------------------------------
    % The absolute direct-path component is modeled with a Rician amplitude.
    % SDR is used to split a cluster's total power between the specular MPC
    % contribution and the DMC contribution.
    p.direct.rician = struct('nu', 0.91, 'sigma', 0.25);
    p.direct.sdr    = struct('mu', 9.90, 'sigma', 2.45);

    % SDR for ordinary MPC clusters, separated by propagation condition.
    p.clusterSDR.los  = struct('mu', 2.19,  'sigma', 3.56);
    p.clusterSDR.nlos = struct('mu', -1.09, 'sigma', 5.42);

    % ---------------------------------------------------------------------
    % 7) Intra-cluster power decay
    % ---------------------------------------------------------------------
    % These parameters control the relative power decay of MPCs within a
    % cluster before and after the delay breakpoint.
    %
    % The 'pre' and 'post' fields correspond to the two regimes in the
    % segmented intra-cluster model.
    p.mpcPowerDecay.los.pre  = struct('zeta', 3.93,   'lambdaTau', 9.31, 'lambdaPhi', 27.65);
    p.mpcPowerDecay.los.post = struct('zeta', -15.26, 'lambdaTau', 2.12, 'lambdaPhi', 3.60);

    p.mpcPowerDecay.nlos.pre  = struct('zeta', -14.25, 'lambdaTau', 8.04, 'lambdaPhi', 7.10);
    p.mpcPowerDecay.nlos.post = struct('zeta', -19.32, 'lambdaTau', 0.06, 'lambdaPhi', 5.96);

    % Direct cluster extra MPC power decay parameters.
    % These are used only for the stochastic MPCs in the direct cluster.
    p.direct.powerDecay.pre  = struct('zeta', 0.79,   'lambdaTau', 24.28, 'lambdaPhi', 35.85);
    p.direct.powerDecay.post = struct('zeta', -30.99, 'lambdaTau', 2.40,  'lambdaPhi', 1.76);

    % ---------------------------------------------------------------------
    % 8) Distance-dependent cluster power models
    % ---------------------------------------------------------------------
    % Cluster power is modeled in dB using a linear function of log-distance,
    % with a distance-dependent standard deviation.
    %
    % The fields mean:
    %   alpha - intercept
    %   beta  - slope in log-distance domain
    %   rho   - baseline standard deviation term
    %   eta   - slope of the standard deviation term
    p.clusterPower.los  = struct('alpha', -22.09, 'beta', -1.35, 'rho', 7.89, 'eta', -0.24);
    p.clusterPower.nlos = struct('alpha', -28.89, 'beta', -2.17, 'rho', 6.97, 'eta', 0.13);

    % Standalone DMC clusters use their own fitted power model.
    p.standaloneDMC.power.los  = struct('alpha', -31.57, 'beta', -1.43, 'rho', 6.57, 'eta', 0.44);
    p.standaloneDMC.power.nlos = struct('alpha', -41.27, 'beta', -1.03, 'rho', 7.39, 'eta', 0.21);

    % ---------------------------------------------------------------------
    % 9) Associated and standalone DMC spread parameters
    % ---------------------------------------------------------------------
    % Associated DMC:
    %   - direct and LoS share the same parameter set
    %   - NLoS has a separate fitted parameter set
    %
    % These parameters control the Gamma-distributed delay/angle decay
    % constants used in the DMC tensor generator.
    p.assocDMC.spread.direct_los.delay = struct('a', 4.87, 'b', 0.27);
    p.assocDMC.spread.direct_los.angle = struct('a', 4.91, 'b', 0.10);
    p.assocDMC.spread.nlos.delay = struct('a', 7.66, 'b', 0.19);
    p.assocDMC.spread.nlos.angle = struct('a', 5.86, 'b', 0.09);

    % Standalone DMC uses its own fitted parameters.
    p.standaloneDMC.spread.los.delay = struct('a', 5.03, 'b', 0.22);
    p.standaloneDMC.spread.los.angle = struct('a', 5.15, 'b', 0.03);
    p.standaloneDMC.spread.nlos.delay = struct('a', 8.31, 'b', 0.16);
    p.standaloneDMC.spread.nlos.angle = struct('a', 5.87, 'b', 0.05);
end