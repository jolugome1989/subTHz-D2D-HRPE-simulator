function [hCluster, info] = d2d_generate_dmc_cluster_tensor( ...
    centroid, spread, clusterPowerLin, grid, cfg, thresholdLin)
% Generate one DMC cluster as a correlated complex random field.
%
% -------------------------------------------------------------------------
% 🔬 MODEL OVERVIEW
% -------------------------------------------------------------------------
% This function generates a diffuse multipath component (DMC) cluster as a
% tensor over delay × Rx × Tx dimensions:
%
%   hCluster ∈ ℂ^(M_f × N_Rx × N_Tx)
%
% where:
%   M_f   → number of delay samples
%   N_Rx  → number of Rx antennas
%   N_Tx  → number of Tx antennas
%
% The generation follows:
%
%   1) Delay-domain correlation (exponential decay model)
%   2) Angular-domain correlation (EADF-based spatial covariance)
%   3) Joint tensor synthesis (Kronecker-type structure)
%   4) Apply delay envelope centered at cluster centroid
%   5) Apply threshold-aware power scaling
%
% -------------------------------------------------------------------------
% ⚠️ THRESHOLD-AWARE POWER SCALING
% -------------------------------------------------------------------------
% Let h0 denote the raw generated tensor.
%
% We seek a scaling factor s such that:
%
%   ∑_{i : |s·h0(i)|² ≥ T} |s·h0(i)|² = P_cluster
%
% where:
%   T = thresholdLin
%   P_cluster = clusterPowerLin
%
% Key properties:
%   - The "visible" bins depend on scaling s
%   - This creates a nonlinear selection problem
%   - Solved via sorting and cumulative energy evaluation
%
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Grid extraction
% -------------------------------------------------------------------------
delayGrid = grid.delay_m(:);      % [M_f × 1]
phiGridR  = grid.phi_R(:);        % Rx angular grid
phiGridT  = grid.phi_T(:);        % Tx angular grid

M_f  = numel(delayGrid);
N_Rx = cfg.N_R;
N_Tx = cfg.N_T;

% -------------------------------------------------------------------------
% 1) Sample DMC spread parameters
% -------------------------------------------------------------------------
% These control exponential decay in delay and angular domains
aDelay  = d2d_sample_gamma(spread.delay.a, spread.delay.b);   % delay spread
lambdaR = d2d_sample_gamma(spread.angle.a, spread.angle.b);   % Rx angular spread
lambdaT = d2d_sample_gamma(spread.angle.a, spread.angle.b);   % Tx angular spread

% -------------------------------------------------------------------------
% 2) Delay-domain correlation
% -------------------------------------------------------------------------
% Build covariance matrix based on exponential decay:
%   R_f(i,j) = exp(-|τ_i - τ_j| / aDelay)
dt  = abs(delayGrid - delayGrid.');
R_f = exp(-dt ./ max(aDelay, eps));

% Ensure Hermitian PSD matrix
R_f = (R_f + R_f')/2 + 1e-9 * eye(M_f);

% Eigen decomposition → generates correlated Gaussian process
[Vf, Df] = eig(R_f);
evals = max(real(diag(Df)), 0);

% Generate uncorrelated complex Gaussian noise
z = (randn(M_f, N_Rx*N_Tx) + 1j*randn(M_f, N_Rx*N_Tx)) / sqrt(2);

% Apply correlation
dmc_delay = reshape(Vf * diag(sqrt(evals)) * z, [M_f, N_Rx, N_Tx]);

% -------------------------------------------------------------------------
% 3) Angular correlation (EADF-based)
% -------------------------------------------------------------------------
% Evaluate array responses on angular grid
psiGrid.phi_R = phiGridR.';
psiGrid.phi_T = phiGridT.';

ArrayPattern = d2d_EvaluateEADFPatternPath(psiGrid, cfg);

delta_phi = max(grid.delta_phi, eps);

% Angular power spectra (Laplacian form)
kR = lambdaR/2 * exp(-lambdaR * abs(d2d_wrap_pi(phiGridR - centroid.aoa_rad))) * delta_phi;
kT = lambdaT/2 * exp(-lambdaT * abs(d2d_wrap_pi(phiGridT - centroid.aod_rad))) * delta_phi;

KR = diag(kR);
KT = diag(kT);

% Spatial covariance matrices
B_R = ArrayPattern.a_R_l * KR * ArrayPattern.a_R_l';
B_T = ArrayPattern.a_T_l * KT * ArrayPattern.a_T_l';

% Ensure Hermitian PSD
B_R = (B_R + B_R')/2 + 1e-9*eye(N_Rx);
B_T = (B_T + B_T')/2 + 1e-9*eye(N_Tx);

% Cholesky decomposition → spatial coloring
L_BR = chol(B_R,'lower');
L_BT = chol(B_T,'lower');

% Combine delay and spatial correlation
h_temp = reshape(reshape(dmc_delay,[M_f*N_Rx,N_Tx]) * L_BT.', [M_f,N_Rx,N_Tx]);

h_temp = permute( ...
    reshape( ...
        reshape(permute(h_temp,[1 3 2]),[M_f*N_Tx,N_Rx]) * L_BR.', ...
        [M_f,N_Tx,N_Rx]), ...
    [1 3 2]);

% -------------------------------------------------------------------------
% 4) Apply delay envelope
% -------------------------------------------------------------------------
% Cluster-centered exponential decay
delayEnv = exp(-abs(delayGrid - centroid.delay_m) / max(aDelay, eps));
h_temp = h_temp .* reshape(delayEnv, [M_f,1,1]);

% -------------------------------------------------------------------------
% 5) Threshold-aware nonlinear scaling
% -------------------------------------------------------------------------
% Flatten power map
p = abs(h_temp(:)).^2;

% Sort descending (largest energy first)
[p_sorted, ~] = sort(p, 'descend');

% Cumulative energy
cumPower = cumsum(p_sorted);

k_valid = 0;
s = 0;

% Search for consistent scaling region
for k = 1:length(p_sorted)

    % Candidate scaling
    s2 = clusterPowerLin / cumPower(k);

    % Check threshold consistency:
    %   s² * p_k ≥ T   AND   s² * p_(k+1) < T
    if s2 * p_sorted(k) >= thresholdLin
        if k == length(p_sorted) || s2 * p_sorted(k+1) < thresholdLin
            k_valid = k;
            s = sqrt(s2);
            break;
        end
    end
end

% Final scaling
if k_valid > 0
    hCluster = h_temp * s;
else
    % Entire cluster falls below threshold → disappears
    hCluster = complex(zeros(size(h_temp)));
end

% -------------------------------------------------------------------------
% 6) Output metadata
% -------------------------------------------------------------------------
info = struct();
info.delayDecay = aDelay;
info.lambdaR = lambdaR;
info.lambdaT = lambdaT;
info.centroid = centroid;
info.clusterPowerLin = clusterPowerLin;

end