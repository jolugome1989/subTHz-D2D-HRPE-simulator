function amp = d2d_sample_rician_amp(nu, sigma)
% Sample a Rician-distributed amplitude.
%
% This function generates a random amplitude following a Rician distribution,
% which is commonly used to model the magnitude of a signal composed of:
%   - a deterministic (LoS) component
%   - a random scattered (NLoS) component
%
% The Rician amplitude is constructed as the magnitude of a complex Gaussian
% random variable with a non-zero mean:
%
%   h = (nu + sigma * N(0,1)) + j * (sigma * N(0,1))
%   amp = |h|
%
% where:
%   - nu    : deterministic (LoS) component magnitude (real part offset)
%   - sigma : standard deviation of the scattered components (both real & imag)
%
% Physically:
%   - large nu / sigma → strong LoS (more deterministic channel)
%   - small nu / sigma → Rayleigh-like fading (pure scattering)
%
% Input:
%   nu    - mean (deterministic) component of the real part
%   sigma - standard deviation of the Gaussian scattering components
%
% Output:
%   amp   - scalar amplitude following a Rician distribution
%
% Example:
%   % Strong LoS case
%   amp = d2d_sample_rician_amp(1.0, 0.2);
%
%   % Rayleigh case (nu = 0)
%   amp = d2d_sample_rician_amp(0, 0.5);

    % ---------------------------------------------------------------------
    % 1) Generate real and imaginary Gaussian components
    % ---------------------------------------------------------------------
    % Real part has a deterministic offset (nu) plus Gaussian noise
    re = nu + sigma * randn;

    % Imaginary part is purely Gaussian
    im = sigma * randn;

    % ---------------------------------------------------------------------
    % 2) Compute magnitude (Rician amplitude)
    % ---------------------------------------------------------------------
    amp = sqrt(re.^2 + im.^2);
end