1. Overview

This repository contains a stochastic HRPE-based channel simulator for sub-THz outdoor D2D communications, modeling both specular (MPC) and diffuse (DMC) propagation components.
It generates high-resolution delay–angle domain channel tensors with physically consistent power and spatial structure.


2. What is needed to run

MATLAB (standard numerical toolboxes are sufficient)
All .m files in this folder added to the MATLAB path
Optional: user-provided EADF matrices (otherwise default isotropic patterns are used)


3. Quick Example
clear all,close all,clc
dbstop if error
tic
d_LOS = 40;
cfg = struct('distance_m', d_LOS, 'scenario', 'LOS');
sim = subthz_d2dhr_simulator(cfg);
Az_axis = 0:10:350;
DDAPS = abs(squeeze(sum(sim.ddir,1))).^2;
[max_power,idx] = max(DDAPS,[],'all');
[idx_r,idx_c] = find(DDAPS == max_power,1);
% PDPs (Max-Dir & Omni)
PDP_temp = abs(sim.ddir).^2;
PDP_md = PDP_temp(:,idx_r,idx_c);
PDP_temp1 =squeeze(max(PDP_temp,[],3));
PDP_omni =squeeze(max(PDP_temp1,[],2)); 
Nf = 1001;
d = (0:Nf-1)*(3e8/1e9); 
figure()
plot(d+d_LOS,10*log10(PDP_md),'b',d+d_LOS,10*log10(PDP_omni),'r'),grid minor, xlabel('Delay (m)'),ylabel('PDP'),legend('MD','Omni')
[xx,yy] = meshgrid(Az_axis);
figure(),surf(xx,yy,10*log10(DDAPS)), shading interp, axis tight, colormap jet, title ('DDAPS'),xlabel('Az_{Tx}'),ylabel('Az_{Rx}')
toc

imagesc(abs(sim.ddir(:,:,1)))
title('Channel magnitude (delay × Rx)')

4. Structure
/simulator
│
├── subthz_d2dhr_simulator.m   % Main simulation script
│
├── MPC generation
│   ├── d2d_generate_mpc_rays.m
│   ├── d2d_make_ray.m
│   ├── d2d_sample_mpc_count.m
│
├── DMC generation
│   ├── d2d_generate_dmc_cluster_tensor.m
│   ├── d2d_make_dmc_cluster.m
│   ├── d2d_make_dmc_cluster_init.m
│
├── Statistical sampling
│   ├── d2d_sample_gamma.m
│   ├── d2d_sample_cluster_power_db.m
│   ├── d2d_sample_cluster_angle.m
│   ├── d2d_sample_rician_amp.m
│   ├── d2d_sample_truncated_cluster_count.m
│
├── Geometry & utilities
│   ├── d2d_segmented_delay_rnd.m
│   ├── d2d_laplace_rnd.m
│   ├── d2d_is_overlapping_standalone.m
│
├── EADF processing
│   ├── d2d_EvaluateEADFPatternPath.m
│   ├── d2d_validate_eadf_dimensions.m
│
├── Configuration
│   ├── d2d_apply_defaults.m
│   ├── d2d_params.m