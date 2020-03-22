function 		d=mult_precond_pcd (g, x, Ax, InsideCGforTN, par,diag_AtA);   
% PCD (Parallel Coordinate Descent) direction  
%
%Call: d=mult_precond_pcd (g, x, Ax, InsideCGforTN, par);   
%
% Input:
%  g - [gradient] vector to be multiplied by preconditioner
%  x, Ax 
%  InsideCGforTN: 1 - when in  process of CG solving Newton system, use just inverse diag of Hessian
%                                  0 - use PCD  (parallel coordinate descent) direction given by 
%  par.func_x      - reference to the user function of x to be optimized
% par.weight_abs_penalty, 
% par.eps_smooth_abs  - arguments of par.weight_abs_penalty, par.eps_smooth_abs
%  par.diag_AtA  -  [approximation of] diag(A'*A) 
%
%Output:
%  d - preconditioned direction

% Michael Zibulevsky 10.07.2008

[f2,g2_x,dummy,diagh2]= par.func_x(x,[],par);


if InsideCGforTN,
	d= g./(diag_AtA+diagh2);   % inverse diag of Hessian

else  % PCD  (parallel coordinate descent) direction
	
	g1=g-g2_x; %grad of quadratic term
	%w=1; % dummy unit weights in least-square objective
	x_s=CoordinateLinesearch_abs_smoothed(x, g1, diag_AtA, par.weight_abs_penalty, par.eps_smooth_abs);
	d=x -x_s;
end

% x_s=CoordinateLinesearch_abs_smoothed(x0,g,w,lambda0,eps);
