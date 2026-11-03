function x = d2d_sample_gamma(a, b)
% Sample from a Gamma distribution.
%
% This function generates random samples from a Gamma distribution with:
%   - shape parameter a
%   - scale parameter b
%
% MATLAB’s built-in function `randg(a)` returns samples from a Gamma
% distribution with:
%   - shape = a
%   - scale = 1 (unit-scale Gamma)
%
% Therefore, multiplying by b converts the samples to the desired
% scale parameter.
%
% Input:
%   a - shape parameter (must be positive)
%   b - scale parameter (must be positive)
%
% Output:
%   x - random sample(s) drawn from Gamma(a, b)
%
% Example:
%   % Generate one sample from Gamma(2, 3)
%   x = d2d_sample_gamma(2, 3);
%
%   % Generate 10 samples (vectorized usage)
%   x = d2d_sample_gamma(2, 3) * ones(10,1);  % or call randg directly

    % Generate unit-scale Gamma samples and apply scaling
    x = randg(a) * b;
end