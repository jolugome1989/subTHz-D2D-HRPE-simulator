function x = d2d_wrap_pi(x)
% Wrap angle(s) to the interval [-pi, pi).
%
% This function maps any real-valued angle (in radians) into the principal
% interval [-pi, pi) using modular arithmetic. It is commonly used for
% azimuth angles, phase differences, and other circular quantities where
% values that differ by integer multiples of 2*pi are equivalent.
%
% Input:
%   x : scalar, vector, or array of angles in radians
%
% Output:
%   x : angle(s) wrapped to [-pi, pi)
%
% Example:
%   -3*pi        -> -pi
%    3*pi/2      -> -pi/2
%    5*pi/4      -> -3*pi/4
%    0           ->  0

    % Shift by pi so that the target interval [-pi, pi) becomes [0, 2*pi),
    % apply modulo operation, then shift back.
    x = mod(x + pi, 2*pi) - pi;
end