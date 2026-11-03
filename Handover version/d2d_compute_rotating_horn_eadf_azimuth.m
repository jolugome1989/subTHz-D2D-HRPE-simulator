function [EADF, meta] = d2d_compute_rotating_horn_eadf_azimuth( ...
    H_meas, phi_meas, phi_steer)
%
% This function computes the azimuth-only EADF (Effective Antenna
% Distribution Function) of a measured horn pattern and then rotates the
% EADF coefficients to emulate a rotating-horn / virtual-array setup.
%
% Key idea:
%   - EADF is the Fourier representation of the measured horn pattern.
%   - Horn steering is emulated in the coefficient domain by multiplying
%     the Fourier coefficients with a phase ramp.
%
% Output convention:
%   - The first dimension corresponds to steering directions, not receive
%     ports.
%   - In other words, N_steer is the number of virtual horn positions.
%
% Supported measurement formats:
%   1) Frequency-independent pattern:
%        H_meas has size [1 x N_ang] or [N_ang x 1]
%      → EADF is returned as [N_steer x N_mode]
%
%   2) Frequency-dependent pattern:
%        H_meas has size [N_f x N_ang]
%      → EADF is returned as [N_steer x N_f x N_mode]
%
% The angular grid can be customized by the user, for example:
%   -180:1:179
%   -90:1:90
%   -60:10:60
%
% If a full 360-degree grid contains a duplicated periodic endpoint
% (for example, -180:1:180 or 0:1:360), the duplicate endpoint is removed.
%
% ========================================================================
% INPUT
% ========================================================================
% H_meas   : measured horn pattern samples
%            Either:
%              [1 x N_ang]          frequency-independent
%              [N_f x N_ang]        frequency-dependent
%
% phi_meas : measured azimuth angles in degrees
%            Examples:
%              -180:1:179
%              -90:1:90
%              -60:10:60
%
% phi_steer : steering directions in degrees
%             Examples:
%               0:10:350
%              -60:10:60
%
% ========================================================================
% OUTPUT
% ========================================================================
% EADF : rotated EADF coefficients
%        If frequency-independent:
%          [N_steer x N_mode]
%        If frequency-dependent:
%          [N_steer x N_f x N_mode]
%
% meta : structure containing metadata
%        Fields include:
%          - type                 : 'frequency-independent' or
%                                    'frequency-dependent'
%          - isFrequencyDependent  : logical flag
%          - N_f                  : number of frequencies / subcarriers
%          - N_ang               : number of cleaned angular samples
%          - N_mode              : number of Fourier modes
%          - N_steer             : number of steering directions
%          - phi_meas_deg        : cleaned measured angular grid
%          - phi_steer_deg       : steering angle grid
%          - mode_index          : centered Fourier mode index vector
%          - rotationConvention  : description of the coefficient rotation
%          - note                : short usage note
%
% ========================================================================


% ------------------------------------------------------------------------
% BASIC INPUT CHECKS
% ------------------------------------------------------------------------
if nargin < 3
    error('Three inputs are required: H_meas, phi_meas, and phi_steer.');
end

if ~isnumeric(H_meas)
    error('H_meas must be numeric.');
end

if ~isvector(phi_meas) || isempty(phi_meas)
    error('phi_meas must be a non-empty vector of measured angles.');
end

if ~isvector(phi_steer) || isempty(phi_steer)
    error('phi_steer must be a non-empty vector of steering angles.');
end


% ------------------------------------------------------------------------
% STANDARDIZE ANGLE INPUTS
% ------------------------------------------------------------------------
% Keep the measured angles in the user-provided coordinate system.
% We only remove a duplicated periodic endpoint if a full 360-degree
% coverage is present.
phi_meas  = phi_meas(:);
phi_steer = phi_steer(:);


% ------------------------------------------------------------------------
% DETECT INPUT FORMAT
% ------------------------------------------------------------------------
% If H_meas is a vector, interpret it as a single frequency-independent
% measured pattern.
if isvector(H_meas)
    H_meas = reshape(H_meas, 1, []);
end

% We support a 2D input [N_f x N_ang] only.
% Each row corresponds to one frequency/subcarrier.
if ndims(H_meas) ~= 2
    error('H_meas must be 2D with size [N_f x N_ang], or a vector of length N_ang.');
end

[N_f, N_ang_in] = size(H_meas);

if N_ang_in ~= numel(phi_meas)
    error('The number of measured angles must match the second dimension of H_meas.');
end

isFrequencyDependent = (N_f > 1);


% ------------------------------------------------------------------------
% SORT THE MEASURED ANGULAR GRID
% ------------------------------------------------------------------------
% Sort the measured angles in ascending order and reorder the measured
% horn pattern to preserve angle-pattern correspondence.
% Note: sorting the angles does not impact Fourier transform
[phi_meas_sorted, sort_idx] = sort(phi_meas, 'ascend');
H_meas = H_meas(:, sort_idx);


% ------------------------------------------------------------------------
% REMOVE DUPLICATED PERIODIC ENDPOINT IF PRESENT
% ------------------------------------------------------------------------
% Example:
%   -180:1:180  -> remove the last 180-degree sample
%   0:1:360     -> remove the last 360-degree sample
%
% This keeps the angular grid periodic without duplicating the same
% direction at both ends.
tol_deg = 1e-12;
if numel(phi_meas_sorted) > 1
    span_deg = phi_meas_sorted(end) - phi_meas_sorted(1);

    if abs(span_deg - 360) < tol_deg
        % Remove the duplicated endpoint.
        phi_meas_sorted(end) = [];
        H_meas(:, end) = [];
    end
end

phi_meas_clean = phi_meas_sorted;
N_ang = numel(phi_meas_clean);


% ------------------------------------------------------------------------
% SET THE FOURIER MODE INDEX VECTOR
% ------------------------------------------------------------------------
% The number of Fourier modes is set equal to the number of cleaned angular
% samples. The mode index is centered around zero.
%
% Examples:
%   N_mode = 360 -> mode_index = -180:179
%   N_mode = 121 -> mode_index = -60:60
N_mode = N_ang;

if mod(N_mode, 2) == 0 % even number
    mode_index = (-N_mode/2) : (N_mode/2 - 1);
else % odd number
    mode_index = (-(N_mode-1)/2) : ((N_mode-1)/2);
end

mode_index = mode_index(:).';   % row vector for matrix construction


% ------------------------------------------------------------------------
% BUILD THE FOURIER BASIS FOR THE MEASURED GRID
% ------------------------------------------------------------------------
% This matrix projects the measured horn pattern into the Fourier domain.
%
% The basis is built from the actual measured angles (in radians), so the
% user can provide custom angular windows such as:
%   -180:1:179
%   -90:1:90
%   -60:10:60
%
% The implementation assumes the measured samples are on a uniform grid.
phi_meas_rad = deg2rad(phi_meas_clean(:));

FourierBasisMeas = 1/sqrt(N_mode) * exp(-1j * (phi_meas_rad * mode_index));


% ------------------------------------------------------------------------
% COMPUTE THE BASE EADF COEFFICIENTS
% ------------------------------------------------------------------------
% The EADF is the Fourier transform of the measured horn pattern.
% Each frequency/subcarrier slice is projected independently.
EADF_base = complex(zeros(N_f, N_mode));

for ff = 1:N_f
    % Measured horn pattern at one frequency/subcarrier:
    % size = [1 x N_ang] after the input normalization above.
    H_slice = H_meas(ff, :);

    % Projection onto the Fourier basis.
    % This produces the EADF coefficients for the selected frequency.
    EADF_base(ff, :) = H_slice * conj(FourierBasisMeas);
end


% ------------------------------------------------------------------------
% BUILD THE STEERING ROTATION OPERATOR
% ------------------------------------------------------------------------
% Horn steering is emulated by rotating the Fourier coefficients.
%
% For steering angle phi_s, the coefficient-domain rotation is:
%   c_rot(m) = c(m) * exp(+j * m * phi_s)
%
% This is the correct Fourier-domain way to emulate the rotated horn
% response in the virtual array setup.
phi_steer_rad = deg2rad(phi_steer);
steerPhase = exp(1j * (mode_index(:) * phi_steer_rad(:).'));   % [N_mode x N_steer]


% ------------------------------------------------------------------------
% APPLY THE STEERING ROTATION
% ------------------------------------------------------------------------
N_steer = numel(phi_steer);

if ~isFrequencyDependent
    % Frequency-independent case:
    % Output size = [N_steer x N_mode]
    EADF = complex(zeros(N_steer, N_mode));

    baseCoeff = EADF_base(1, :).';   % [N_mode x 1]

    % Apply the steering phase to each steering direction.
    % Result per steering angle: [1 x N_mode]
    EADF = (baseCoeff .* steerPhase).';

else
    % Frequency-dependent case:
    % Output size = [N_steer x N_f x N_mode]
    EADF = complex(zeros(N_steer, N_f, N_mode));

    for ff = 1:N_f
        baseCoeff = EADF_base(ff, :).';   % [N_mode x 1]

        % Rotate coefficients for all steering directions.
        % The transpose gives [N_steer x N_mode].
        rotated = (baseCoeff .* steerPhase).';

        % Store as [N_steer x 1 x N_mode] inside the 3D output.
        EADF(:, ff, :) = reshape(rotated, N_steer, 1, N_mode);
    end
end


% ------------------------------------------------------------------------
% METADATA
% ------------------------------------------------------------------------
meta = struct();
if isFrequencyDependent
    meta.type = 'frequency-dependent';
else
    meta.type = 'frequency-independent';
end

meta.isFrequencyDependent = isFrequencyDependent;
meta.N_f = N_f;
meta.N_ang = N_ang;
meta.N_mode = N_mode;
meta.N_steer = N_steer;
meta.phi_meas_deg = phi_meas_clean(:).';
meta.phi_steer_deg = phi_steer(:).';
meta.mode_index = mode_index;
meta.rotationConvention = 'c_rot(m) = c(m) * exp(+j * m * phi_steer_rad)';
meta.note = ['EADF is the Fourier transform of the measured horn ' ...
             'pattern. The output EADF is already rotated for each ' ...
             'steering direction.'];

end