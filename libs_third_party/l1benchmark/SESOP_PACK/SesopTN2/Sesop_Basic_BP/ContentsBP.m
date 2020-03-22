%                     Example of use SesopTn:  Basis Pursuit
%
%            Find sparse solution of underdetermined linear system  Ax=b
%            with noisy right-hand side  y=b+noise, 
%            using smoothed 1-norm penalty: 
%
%                    min_x  ||Ax-y||^2 +  mu*sum(abs_smooth(x))
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Files:
%
%   BasisPursuitSesop                -  simple script
%   BasisPursuitSesop_paper_fig.m    -  advanced script
%
%   mult_diag_precond    -  Diagonal preconditioner: d = g./diag(Hessian)


% Michael Zibulevsky 10.07.2008;   06.08.2008; 03.02.2009
%
% Copyright (c) 2008. Free for academic use. No warranty 
