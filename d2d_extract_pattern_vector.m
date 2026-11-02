function a = d2d_extract_pattern_vector(ArrayPattern, side, N)
% Extract a single-path array response vector from EADF output.
%
% This function converts the output of EvaluateEADFPatternPath into a
% clean column vector representing the array response for a single path.
%
% Since the EADF function can return data in different shapes depending on
% configuration (e.g., frequency-dependent, multi-path, or singleton),
% this function:
%   - selects the correct field (Rx or Tx)
%   - removes unnecessary dimensions
%   - ensures the result is a column vector of size [N x 1]
%
% Input:
%   ArrayPattern - struct returned by EvaluateEADFPatternPath, possibly
%                  containing fields:
%                  a_R_l, a_T_l, a_R, a_T, etc.
%
%   side         - 'R' for receiver or 'T' for transmitter
%
%   N            - expected number of antenna elements
%
% Output:
%   a            - column vector [N x 1] representing array response
%
% Example:
%   aR = d2d_extract_pattern_vector(ArrayPattern, 'R', 36);
%   aT = d2d_extract_pattern_vector(ArrayPattern, 'T', 36);

    % ---------------------------------------------------------------------
    % 1) Select appropriate field based on side
    % ---------------------------------------------------------------------
    switch upper(side)

        case 'R'  % Receiver side
            if isfield(ArrayPattern, 'a_R_l')
                tmp = ArrayPattern.a_R_l;
            elseif isfield(ArrayPattern, 'a_R')
                tmp = ArrayPattern.a_R;
            else
                % Fallback: isotropic response
                tmp = ones(N, 1);
            end

        case 'T'  % Transmitter side
            if isfield(ArrayPattern, 'a_T_l')
                tmp = ArrayPattern.a_T_l;
            elseif isfield(ArrayPattern, 'a_T')
                tmp = ArrayPattern.a_T;
            else
                tmp = ones(N, 1);
            end

        otherwise
            error('Unknown side "%s".', side);
    end

    % ---------------------------------------------------------------------
    % 2) Remove extra dimensions (robust handling)
    % ---------------------------------------------------------------------
    % Possible shapes:
    %   [N x M_f x L]  -> take first frequency & path
    %   [N x L]        -> take first path
    %   [N x 1]        -> already correct
    if ndims(tmp) == 3
        tmp = tmp(:,1,1);        % take first frequency and path
    elseif ismatrix(tmp)
        if size(tmp, 2) > 1
            tmp = tmp(:,1);      % take first column
        else
            tmp = tmp(:);        % ensure column vector
        end
    else
        tmp = tmp(:);
    end

    % ---------------------------------------------------------------------
    % 3) Final safety check
    % ---------------------------------------------------------------------
    % If empty, fallback to isotropic response
    if isempty(tmp)
        a = ones(N, 1);
    else
        a = tmp(:);  % enforce column vector format
    end
end