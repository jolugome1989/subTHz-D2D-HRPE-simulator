function [ArrayPattern] = d2d_EvaluateEADFPatternPath(psi, ParaSounder)
% Evaluate array response using EADF (Effective Aperture Distribution Function).
%
% -------------------------------------------------------------------------
%   REFERENCES
% -------------------------------------------------------------------------
% [1] M. Landmann,
%     "Efficient Antenna Description for MIMO Channel Modelling and Estimation,"
%     European Conference on Wireless Technology (ECWT), 2004.
%
% -------------------------------------------------------------------------
%   THEORY: EADF REPRESENTATION
% -------------------------------------------------------------------------
% The antenna response is represented via Fourier expansion:
%
%   a(φ, θ) = Σ_k Σ_l c(k,l) · exp(-j k φ) · exp(-j l θ)
%
% Azimuth-only:
%   a(φ) = Σ_k c(k) · exp(-j k φ)
%
% where:
%   c(k,l) → EADF coefficients
%
% -------------------------------------------------------------------------
%   NORMALIZATION FACTOR
% -------------------------------------------------------------------------
% The normalization factor ensures consistency between DFT and inverse DFT.
%
% In this implementation:
%
%   f(φ) = (1 / sqrt(N_modes)) * exp(-j k φ)
%
% where:
%   N_modes = number of Fourier coefficients = (2K+1)
%
% Practical interpretation:
%   If no compression is applied:
%       1° sampling → 360 samples → NormFactor = 360
%       5° sampling → 72 samples → NormFactor = 72
%
% IMPORTANT:
%   NormFactor = number of retained Fourier modes (NOT always raw samples)
%
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% 1) Normalization factors
% -------------------------------------------------------------------------
if isfield(ParaSounder,'NormFactorPhi_RX')==0
    NormFactorPhi_RX = size(ParaSounder.EADF_RX, ndims(ParaSounder.EADF_RX));
else
    NormFactorPhi_RX = ParaSounder.NormFactorPhi_RX;
end

if isfield(ParaSounder,'NormFactorPhi_TX')==0
    NormFactorPhi_TX = size(ParaSounder.EADF_TX, ndims(ParaSounder.EADF_TX));
else
    NormFactorPhi_TX = ParaSounder.NormFactorPhi_TX;
end

% Default sampling mode
if isfield(ParaSounder,'sampling_mode')==0
    ParaSounder.sampling_mode = 'uniform';
end

ArrayPattern.temp = [];

c = 299792458; % speed of light
M_f = ParaSounder.M_f;
N_R = ParaSounder.N_R;
N_T = ParaSounder.N_T;

% -------------------------------------------------------------------------
% 2) Main processing
% -------------------------------------------------------------------------
switch ParaSounder.sampling_mode

    % =====================================================================
    % CASE: UNIFORM SAMPLING
    % =====================================================================
    case 'uniform'

        switch ParaSounder.ArrayFrequencyDependence

            % =============================================================
            % CASE: Frequency-dependent EADF
            % =============================================================
            case 'FrequencyDependent'

                % ---------------- Receiver side ----------------
                if N_R > 1
                    switch ParaSounder.estimationmode

                        case 'azimuth'

                            L = length(psi.phi_R);

                            ArrayPattern.a_R_l = zeros(N_R, M_f, L);
                            ArrayPattern.a_R_dot_l_phi = zeros(N_R, M_f, L);

                            % Extract number of Fourier modes
                            [~,~,N_EADF_phi] = size(ParaSounder.EADF_RX);
                            N_EADF_phi = (N_EADF_phi - 1)/2;

                            for m = 1:N_R

                                % Extract EADF coefficients for antenna m
                                temp = squeeze(ParaSounder.EADF_RX(m,:,:));

                                for l = 1:L

                                    % Fourier basis exp(-j*k*φ)
                                    f_extrapol_phi = 1/sqrt(NormFactorPhi_RX) * ...
                                        exp(-1j * psi.phi_R(l) * (-N_EADF_phi:N_EADF_phi));

                                    % Derivative basis
                                    f_extrapol_dot_phi = 1/sqrt(NormFactorPhi_RX) * ...
                                        (-1j*(-N_EADF_phi:N_EADF_phi)) .* ...
                                        exp(-1j * psi.phi_R(l) * (-N_EADF_phi:N_EADF_phi));

                                    % Fourier reconstruction
                                    ArrayPattern.a_R_l(m,:,l) = temp * f_extrapol_phi.';
                                    ArrayPattern.a_R_dot_l_phi(m,:,l) = temp * f_extrapol_dot_phi.';
                                end
                            end
                    end
                end

                % ---------------- Transmitter side ----------------
                if N_T > 1
                    switch ParaSounder.estimationmode

                        case 'azimuth'

                            L = length(psi.phi_T);

                            ArrayPattern.a_T_l = zeros(N_T, M_f, L);
                            ArrayPattern.a_T_dot_l_phi = zeros(N_T, M_f, L);

                            [~,~,N_EADF_phi] = size(ParaSounder.EADF_TX);
                            N_EADF_phi = (N_EADF_phi - 1)/2;

                            for m = 1:N_T
                                temp = squeeze(ParaSounder.EADF_TX(m,:,:));

                                for l = 1:L

                                    f_extrapol_phi = 1/sqrt(NormFactorPhi_TX) * ...
                                        exp(-1j * psi.phi_T(l) * (-N_EADF_phi:N_EADF_phi));

                                    f_extrapol_dot_phi = 1/sqrt(NormFactorPhi_TX) * ...
                                        (-1j*(-N_EADF_phi:N_EADF_phi)) .* ...
                                        exp(-1j * psi.phi_T(l) * (-N_EADF_phi:N_EADF_phi));

                                    ArrayPattern.a_T_l(m,:,l) = temp * f_extrapol_phi.';
                                    ArrayPattern.a_T_dot_l_phi(m,:,l) = temp * f_extrapol_dot_phi.';
                                end
                            end
                    end
                end

            % =============================================================
            % CASE: Frequency-independent EADF
            % =============================================================
            case 'FrequencyIndependent'

                if N_R > 1

                    L = length(psi.phi_R);

                    ArrayPattern.a_R_l = zeros(N_R, L);

                    % Extract Fourier modes
                    [~,N_EADF_phi] = size(ParaSounder.EADF_RX);
                    N_EADF_phi = (N_EADF_phi - 1)/2;

                    for m = 1:N_R
                        temp = squeeze(ParaSounder.EADF_RX(m,:));

                        for l = 1:L

                            % Fourier basis
                            f_extrapol_phi = 1/sqrt(NormFactorPhi_RX) * ...
                                exp(-1j * psi.phi_R(l) * (-N_EADF_phi:N_EADF_phi));

                            % Reconstruction
                            ArrayPattern.a_R_l(m,l) = f_extrapol_phi * temp.';
                        end
                    end
                end

                % -------- Transmitter side --------
                if N_T > 1

                    L = length(psi.phi_T);

                    ArrayPattern.a_T_l = zeros(N_T, L);

                    [~,N_EADF_phi] = size(ParaSounder.EADF_TX);
                    N_EADF_phi = (N_EADF_phi - 1)/2;

                    for m = 1:N_T
                        temp = squeeze(ParaSounder.EADF_TX(m,:));

                        for l = 1:L

                            f_extrapol_phi = 1/sqrt(NormFactorPhi_TX) * ...
                                exp(-1j * psi.phi_T(l) * (-N_EADF_phi:N_EADF_phi));

                            ArrayPattern.a_T_l(m,l) = f_extrapol_phi * temp.';
                        end
                    end
                end
        end
end

% -------------------------------------------------------------------------
% SUMMARY
% -------------------------------------------------------------------------
% This function converts EADF coefficients into physical array responses:
%
%   a(φ) = Σ_k c(k) · exp(-j k φ)
%
% The normalization ensures:
%   - energy consistency
%   - correct reconstruction of measured patterns
%
% This function is the core link between:
%   antenna modeling (EADF)
%   and channel simulation (MPC/DMC)
% -------------------------------------------------------------------------

end