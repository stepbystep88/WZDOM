% Basis Pursuit via PCD-Sesop (parallel coordinate descent) using SesopTN optimization tool
%
%          min_c  ||Ac-y||^2 +  mu*sum(smooth_abs(c))
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% SCRIPTS:
%
%   BasisPursuitSesopPCD         - Basis Pursuit using SesopTn and PCD (parallel coordinate descent)
%   BasisPursuitSesopPCDconcave  - "Concave norm" Basis Pursuit using SesopTn and PCD
%
% FUNCTIONS:
%
%   mult_diag_precond       - Diagonal preconditioner: d = g./diag(Hessian)
%   mult_precond_pcd        - PCD (Parallel Coordinate Descent) direction  
%   CoordinateLinesearch_abs_smoothed - Minimize   0.5w*(x_s - x0)^2 +g*x_s+ lambda0*smooth_abs(x_s,eps)
%
%
%  See also: M. Elad, B. Matalon, and M. Zibulevsky, "Coordinate and Subspace Optimization Methods for 
%  Linear Least Squares with Non-Quadratic Regularization", 
%  Applied and Computational Harmonic Analysis, Vol. 23, pp. 346-367, November 2007


%  Michael  Zibulevsky  04.08.2008 04.02.2009
%
% Copyright (c) 2008-2009. All rights reserved. No warranty. Free for academic use

