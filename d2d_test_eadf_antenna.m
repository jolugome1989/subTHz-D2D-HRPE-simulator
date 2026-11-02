clear all,close all,clc
dbstop if error
tic
load("D:\WiDES\SuperSolution\Zihang_Handover\sub-THz\HR_modeling\3_simulator\EADF_az_pattern_2025Jul.mat");
phi_meas = (0:1:358).';
phi_steer = (0:10:350).';
% figure(),mesh(20*log10(abs(EADF_az1))),grid minor, shading interp, colormap jet, axis tight
% figure(),mesh(20*log10(abs(EADF_az2))),grid minor, shading interp, colormap jet, axis tight
G_az1 = ifft(EADF_az1,[],2);
G_az2 = ifft(EADF_az2,[],2);
G_az1_ave = squeeze(mean(G_az1,1));
G_az2_ave = squeeze(mean(G_az2,1));
figure(),plot(phi_meas, 20*log10(abs(G_az1_ave))),grid minor,xlabel('Az'),ylabel('El')
figure(),plot(phi_meas, 20*log10(abs(G_az2_ave))),grid minor,xlabel('Az'),ylabel('El')


% Convert measured EADF data to evaluator-ready rotating-horn patterns.
[EADF1, meta1] = d2d_compute_rotating_horn_eadf_azimuth(G_az1_ave, phi_meas, phi_steer);
[EADF2, meta2] = d2d_compute_rotating_horn_eadf_azimuth(G_az2_ave, phi_meas, phi_steer);
toc