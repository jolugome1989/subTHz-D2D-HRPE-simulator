function x = d2d_wrap_2pi(x)
% Wrap angle(s) to the interval [0, 2*pi).
%
% This function maps any real-valued angle (in radians) into the principal
% interval [0, 2*pi) using modular arithmetic. It is commonly used for
% azimuth angles, phases, and other circular quantities where values that
% differ by integer multiples of 2*pi are equivalent.
%
% Input:
%   x : scalar, vector, or array of angles in radians
%
% Output:
%   x : angle(s) wrapped to [0, 2*pi)
%
% Example:
%   -2*pi/3   ->  4*pi/3
%    0        ->  0
%    7*pi/2   ->  3*pi/2

    % Reduce the input modulo 2*pi.
    x = mod(x, 2*pi);
end