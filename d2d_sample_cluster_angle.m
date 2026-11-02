function ang = d2d_sample_cluster_angle(model)
% Sample a cluster centroid angle.
%
% This function generates a random azimuth angle for a cluster centroid
% based on the specified angular distribution model. Two types of models
% are supported:
%
%   1) 'bimodal_laplace':
%      - A mixture of two Laplace (double-exponential) distributions.
%      - Each component has its own center (mu) and spread (lambda).
%      - The mixture weight vector 'w' determines which component is used.
%      - This models directional clustering with two dominant lobes.
%
%   2) 'uniform':
%      - Uniform distribution over [0, 2*pi).
%      - Used to represent isotropic or highly scattered environments.
%
% The output angle is always wrapped to [0, 2*pi).
%
% Input:
%   model - struct describing the angular distribution:
%           For 'bimodal_laplace':
%               model.mode    = 'bimodal_laplace'
%               model.w       = [w1, w2] mixture weights (sum to 1)
%               model.mu      = [mu1, mu2] component centers (radians)
%               model.lambda  = [b1, b2] Laplace scale parameters
%
%           For 'uniform':
%               model.mode    = 'uniform'
%
% Output:
%   ang   - sampled angle in radians, wrapped to [0, 2*pi)
%
% Example:
%   model.mode = 'bimodal_laplace';
%   model.w = [0.5, 0.5];
%   model.mu = [0, pi];
%   model.lambda = [0.2, 0.8];
%   ang = d2d_sample_cluster_angle(model);

    switch lower(model.mode)

        case 'bimodal_laplace'
            % -------------------------------------------------------------
            % 1) Select one of the two Laplace components
            % -------------------------------------------------------------
            % Draw a Bernoulli trial using the mixture weight
            if rand < model.w(1)
                mu = model.mu(1);
                b  = model.lambda(1);
            else
                mu = model.mu(2);
                b  = model.lambda(2);
            end

            % -------------------------------------------------------------
            % 2) Sample from the selected Laplace distribution
            % -------------------------------------------------------------
            % Truncated to [-pi, pi] to avoid extreme angular outliers
            ang = d2d_wrap_2pi(d2d_laplace_rnd(mu, b, 1, pi));

        case 'uniform'
            % -------------------------------------------------------------
            % Uniform distribution over the full angular range
            % -------------------------------------------------------------
            ang = 2 * pi * rand;

        otherwise
            % -------------------------------------------------------------
            % Error handling for unsupported models
            % -------------------------------------------------------------
            error('Unknown angle model "%s".', model.mode);
    end
end