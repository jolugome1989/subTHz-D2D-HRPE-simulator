function d2d_validate_eadf_dimensions(cfg)
% Validate EADF matrix dimensions based on evaluation mode.
%
% This function checks consistency between the provided EADF matrices and
% the simulator configuration, taking into account:
%   - frequency dependence
%   - sampling mode
%   - estimation mode (azimuth vs full 2D)
%
% The expected dimensionality of EADF depends on:
%   cfg.ArrayFrequencyDependence
%   cfg.sampling_mode
%   cfg.estimationmode
%
% Input:
%   cfg - configuration struct containing:
%         EADF_RX, EADF_TX
%         N_Rx, N_Tx
%         M_f
%         sampling_mode
%         ArrayFrequencyDependence
%         estimationmode
%
% Output:
%   (none) – throws error if mismatch detected

    % ---------------------------------------------------------------------
    % 1) Validate RX EADF
    % ---------------------------------------------------------------------
    if ~isempty(cfg.EADF_RX)

        sz = size(cfg.EADF_RX);

        switch cfg.ArrayFrequencyDependence

            % =============================================================
            % CASE 1: Frequency Independent
            % =============================================================
            case 'FrequencyIndependent'

                if strcmp(cfg.estimationmode, 'azimuth')
                    % Expected: [N_Rx × (2K+1)]
                    if ndims(cfg.EADF_RX) ~= 2
                        error('EADF_RX must be 2D for FrequencyIndependent + azimuth.');
                    end
                    if sz(1) ~= cfg.N_Rx
                        error('EADF_RX dimension mismatch: first dim must be N_Rx.');
                    end

                elseif strcmp(cfg.estimationmode, 'both')
                    % Expected: [N_Rx × (2K_phi+1) × (2K_theta+1)]
                    if ndims(cfg.EADF_RX) ~= 3
                        error('EADF_RX must be 3D for FrequencyIndependent + both.');
                    end
                    if sz(1) ~= cfg.N_Rx
                        error('EADF_RX dimension mismatch: first dim must be N_Rx.');
                    end
                end

            % =============================================================
            % CASE 2: Frequency Dependent
            % =============================================================
            case 'FrequencyDependent'

                if strcmp(cfg.estimationmode, 'azimuth')
                    % Expected: [N_Rx × M_f × (2K+1)]
                    if ndims(cfg.EADF_RX) ~= 3
                        error('EADF_RX must be 3D for FrequencyDependent + azimuth.');
                    end
                    if sz(1) ~= cfg.N_Rx || sz(2) ~= cfg.M_f
                        error('EADF_RX dimension mismatch: [N_Rx × M_f × ...]');
                    end

                elseif strcmp(cfg.estimationmode, 'both')
                    % Expected: [N_Rx × M_f × (2K_phi+1) × (2K_theta+1)]
                    if ndims(cfg.EADF_RX) ~= 4
                        error('EADF_RX must be 4D for FrequencyDependent + both.');
                    end
                    if sz(1) ~= cfg.N_Rx || sz(2) ~= cfg.M_f
                        error('EADF_RX dimension mismatch: [N_Rx × M_f × ...]');
                    end
                end
                
            otherwise
                error('Unknown ArrayFrequencyDependence mode.');
        end
    end

    % ---------------------------------------------------------------------
    % 2) Validate TX EADF (same logic)
    % ---------------------------------------------------------------------
    if ~isempty(cfg.EADF_TX)

        sz = size(cfg.EADF_TX);

        switch cfg.ArrayFrequencyDependence

            case 'FrequencyIndependent'

                if strcmp(cfg.estimationmode, 'azimuth')
                    if ndims(cfg.EADF_TX) ~= 2
                        error('EADF_TX must be 2D for FrequencyIndependent + azimuth.');
                    end
                    if sz(1) ~= cfg.N_Tx
                        error('EADF_TX dimension mismatch: first dim must be N_Tx.');
                    end

                elseif strcmp(cfg.estimationmode, 'both')
                    if ndims(cfg.EADF_TX) ~= 3
                        error('EADF_TX must be 3D for FrequencyIndependent + both.');
                    end
                    if sz(1) ~= cfg.N_Tx
                        error('EADF_TX dimension mismatch: first dim must be N_Tx.');
                    end
                end

            case 'FrequencyDependent'

                if strcmp(cfg.estimationmode, 'azimuth')
                    if ndims(cfg.EADF_TX) ~= 3
                        error('EADF_TX must be 3D for FrequencyDependent + azimuth.');
                    end
                    if sz(1) ~= cfg.N_Tx || sz(2) ~= cfg.M_f
                        error('EADF_TX dimension mismatch.');
                    end

                elseif strcmp(cfg.estimationmode, 'both')
                    if ndims(cfg.EADF_TX) ~= 4
                        error('EADF_TX must be 4D for FrequencyDependent + both.');
                    end
                    if sz(1) ~= cfg.N_Tx || sz(2) ~= cfg.M_f
                        error('EADF_TX dimension mismatch.');
                    end
                end

            otherwise
                error('Unknown ArrayFrequencyDependence mode.');
        end
    end
end