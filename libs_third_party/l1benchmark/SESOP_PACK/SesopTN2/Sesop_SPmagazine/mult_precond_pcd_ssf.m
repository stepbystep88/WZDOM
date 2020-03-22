function 		d=mult_precond_pcd_ssf(g, x, Ax, InsideCGforTN, par);   
% PCD or SSF direction (Parallel Coordinate Descent or Separable Surrogate Function) 
%
%Call: d=mult_precond_pcd_ssf (g, x, Ax, InsideCGforTN, par);   
%
% Input:
%  g - [gradient] vector to be multiplied by preconditioner
%  x, Ax 
%  InsideCGforTN: 1 - when in  process of CG solving Newton system, use just inverse diag of Hessian
%                 0 - use PCD/SSF  (parallel coordinate descent) direction given by...
%
%  par.FlagSSF:   1 -SSF, 0 -PCD
%  par.func_x   - reference to the user function of x to be optimized
%  par.weight_abs_penalty, 
%  par.eps_smooth_abs  - arguments of par.weight_abs_penalty, par.eps_smooth_abs
%  par.diag_AtA  -  [approximation of] diag(A'*A) 
%
%Output:
%  d - preconditioned direction

% Michael Zibulevsky 10.07.2008; 25.02.2010

persistent diag_of_precond

[f2,g2_x,dummy,diagh2]= par.func_x(x,[],par);


if InsideCGforTN,
   if par.FlagPCDSSF_inside_TN && par.eps_smooth_abs>1e-20,
      d= g./diag_of_precond;
   else
      d= g./(par.diag_AtA+diagh2);   % inverse diag of Hessian
   end
else  % PCD  (parallel coordinate descent) or SSF (Separable Surrogate Function) direction 
	
	g1=g-g2_x; %grad of quadratic term
	%w=1; % dummy unit weights in least-square objective
   if par.FlagSSF, ddd=((1.01*par.max_sing_value_A)^2)*ones(numel(g),1); % Quadratic term of Separable Surrogate Function
   else            ddd=par.diag_AtA;   % PCD
   end
	x_s=CoordinateLinesearch_abs_smoothed(x, g1, ddd, par.weight_abs_penalty, par.eps_smooth_abs);
	d=x -x_s;
   if par.eps_smooth_abs>1e-20, diag_of_precond= (g./d)+1e-3;end
end

% x_s=CoordinateLinesearch_abs_smoothed(x0,g,w,lambda0,eps);
