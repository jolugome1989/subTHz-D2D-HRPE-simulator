function field = d2d_path_field_from_ray(ray, cfg)
% Compute the spatial channel contribution of a single ray.
%
% This function evaluates the array response of one MPC (or direct path)
% using the EADF-based antenna model, and constructs the corresponding
% MIMO channel matrix contribution.
%
% The resulting field represents a rank-1 channel component:
%   H_ray = sqrt(P) * exp(j*phi) * (a_R * a_T^H)
%
% where:
%   - a_R : receive array response vector
%   - a_T : transmit array response vector
%   - P   : ray power
%   - phi : random phase
%
% Input:
%   ray - struct containing ray parameters:
%         aoa_rad   : angle of arrival (radians)
%         aod_rad   : angle of departure (radians)
%         power_lin : linear-scale power
%         phase_rad : random phase (radians)
%
%   cfg - configuration struct containing array/EADF settings:
%         N_R       : number of Rx elements
%         N_T       : number of Tx elements
%         (and fields required by d2d_EvaluateEADFPatternPath)
%
% Output:
%   field - complex matrix of size [N_R x N_T]
%           representing the spatial response of this ray
%
% Example:
%   field = d2d_path_field_from_ray(ray, cfg);

    % ---------------------------------------------------------------------
    % 1) Construct input for EADF evaluation
    % ---------------------------------------------------------------------
    % psi contains the angular parameters at which the array response is
    % evaluated (single path case)
    psi = struct();
    psi.phi_R = ray.aoa_rad;
    psi.phi_T = ray.aod_rad;

    % ---------------------------------------------------------------------
    % 2) Evaluate array patterns using EADF
    % ---------------------------------------------------------------------
    ArrayPattern = d2d_EvaluateEADFPatternPath(psi, cfg);

    % ---------------------------------------------------------------------
    % 3) Extract receive and transmit array response vectors
    % ---------------------------------------------------------------------
    aR = d2d_extract_pattern_vector(ArrayPattern, 'R', cfg.N_R);
    aT = d2d_extract_pattern_vector(ArrayPattern, 'T', cfg.N_T);

    % ---------------------------------------------------------------------
    % 4) Construct the rank-1 channel contribution
    % ---------------------------------------------------------------------
    % Outer product a_R * a_T^H gives spatial response
    field = sqrt(ray.power_lin) * exp(1j * ray.phase_rad) * (aR * aT.');
end