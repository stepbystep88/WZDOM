% In order to approximately invert Radon transfotm, we compute
% reconstruction functional "w", which has minimal L2-norm given permitted bluring 
% window "wind". Penalty "lam" compromises between ||w|| and the
% influence of pixels located out of the window.
%
%---------------------- By Michael Zibulevsky, 28.03.2006; version 14.04.2006 ------------------------------%
%
% Contents.m  -  this short readme file 
% compute_w  - the main script
% myradon, myadjradon - Radon transform and its adjoint. These are "slow" functions, fast versions should come soon
% CGmz  - conjugate gradient solver of a symmetric positive semidefinite linear system of equations Hw=b 
% fw        - computing Hw to be solved with CG
% report  - plotting results of the current iteration of CG
%
%
% I also put here an example of ordinary image reconstruction by min ||Rx-y||^2 using CG:
%
% tomogr_recon - main script
% RtR                 - computing R'Rx
% report_recon  - reporting iterations CG of the reconstruction

