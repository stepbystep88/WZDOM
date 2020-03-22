% Compressive Sensing, Image Deblurring and Tomography via L1-L2 optimization 
% with various methods:
%
%  - PCD(Parallel Coordinate Descent) and SSF(Separable Surrogate Functions), 
%    combined with SESOP or CG 
%  - FISTA,  by Beck & Teboulle
%  - L-BFGS, implemented by Mark Schmidt
%  - L1-L2 Interior Point Method,  by Boyd et al.
%
% See the paper 
% "L1-L2 Optimization in Signal and Image Processing: Iterative Shrinkage and Beyond"
% by Michael Zibulevsky and Michael Elad, IEEE Signal Processing Magazine, to appear
% 
% Michael Zibulevsky, 26.02.2010
%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Scripts
%   BasisPursuitComprSens             - L1-L2 Compressive Sensing 
%   BasisPursuitDeblur                - L1-L2 Image Deblurring 
%   BasisPursuitTomogr                - L1-L2 Tomography 
%
%   

% Auxiliary scripts
%   create_test_image                 - Create/read test image X00 and resize it to [m0,n0]
%   init_BP                           - init variables
%   optimize_all                      - run various optimization methods:PCD and SSF combined with SESOP;
%                                         FISTA; L-BFGS, Interior Point
%   OptimizeBP                        - OptimizeBP
%   plot_all                          - plot_all.m
%   show_all_results_FINAL3           - Show results  of all optimization methods

% Functions
%   blurring                          - Convolve signal/image with blur kernel 
%   blurring_adj                      - Adjoint convolution
%   ComprSensProj                     - Compute sub-sampled 2d FFT (indeces in par.ComprSens.ind)
%   ComprSensProj_adj                 - Its adjoint
%   CoordinateLinesearch_abs_smoothed - Minimize in x_s        0.5w*(x_s - x0)^2 +g*x_s+ lambda0*phi(x_s,eps)
%   deblur_report_func                - plot current optimization results 
%   mult_diag_precond                 - Diagonal preconditioner: d = g./diag(Hessian)
%   mult_precond_pcd_ssf              - PCD or SSF direction (Parallel Coordinate Descent or Separable Surrogate Function) 
%   projector_synthesis               - Apply synthesis operator (e.g. wavelets) and then projector (e.g. Blur, Radon)
%   analysis_projector_adj            - Adjoint to projector_synthesis
