function ok = d2d_is_overlapping_standalone(dmcCentroid, occupiedCentroids, cfg)
% Check whether a standalone DMC centroid overlaps with existing clusters.
%
% This function determines if a newly generated standalone DMC cluster
% centroid is too close to any previously established cluster centers.
%
% The set "occupiedCentroids" includes:
%   - MPC cluster centers (implicitly represented by associated DMC centers,
%     since associated DMC shares the same centroid as its parent MPC cluster)
%   - previously accepted standalone DMC cluster centers
%
% Therefore, this function enforces separation between:
%   - standalone DMC vs MPC clusters (via associated DMC centroids)
%   - standalone DMC vs previously generated standalone DMC clusters
%
% The overlap criterion is defined jointly in delay and angular domains:
%   - delay difference below cfg.separationTau_m
%   - AoA difference below cfg.separationAngle_deg
%   - AoD difference below cfg.separationAngle_deg
%
% If the new standalone DMC centroid is too close to ANY existing centroid
% (MPC or DMC), it is considered overlapping and should be rejected and
% regenerated in the main simulator loop.
%
% Input:
%   dmcCentroid        - struct with fields:
%                        delay_m, aoa_rad, aod_rad
%   occupiedCentroids  - array of structs representing existing cluster
%                        centers (MPC centers + accepted standalone DMC centers)
%   cfg                - configuration struct containing separation thresholds:
%                        separationTau_m, separationAngle_deg
%
% Output:
%   ok : boolean
%        true  -> overlap detected (reject and regenerate)
%        false -> no overlap (accept this standalone DMC cluster)
%
% Example:
%   if d2d_is_overlapping_standalone(newCentroid, occupiedCentroids, cfg)
%       % regenerate this standalone DMC cluster
%   else
%       % accept and append to occupiedCentroids
%   end

    ok = false;

    for k = 1:numel(occupiedCentroids)
        % Compute delay difference (in meters)
        dt = abs(dmcCentroid.delay_m - occupiedCentroids(k).delay_m);

        % Compute angular differences (wrapped to [-pi, pi))
        da = abs(d2d_wrap_pi(dmcCentroid.aoa_rad - occupiedCentroids(k).aoa_rad));
        dd = abs(d2d_wrap_pi(dmcCentroid.aod_rad - occupiedCentroids(k).aod_rad));

        % Check joint overlap condition in delay and angular domains
        if dt < cfg.separationTau_m && ...
           rad2deg(da) < cfg.separationAngle_deg && ...
           rad2deg(dd) < cfg.separationAngle_deg

            % Overlap detected: reject this candidate
            ok = true;
            return;
        end
    end
end