function 		d=mult_diag_precond (g, x, Ax, InsideCGforTN, par,diag_AtA);   
% Diagonal preconditioner: d = g./diag(Hessian)
%
%Call: d=mult_diag_precond (g, x, par)
%
% Input:
%  g - [gradient] vector to be multiplied by preconditioner
%  x, par - arguments used  in  Hessian diagonal calculation by user function par.func_x;
%  par.func_x      - reference to the user function of x to be optimized
%  diag_AtA  -  [approximation of] diag(A'*A) 
%
%Output:
%  d - preconditioned direction

% Michael Zibulevsky 10.07.2008

[f2,g2_x,dummy,diagh2]= par.func_x(x,[],par);
d= g./(diag_AtA+diagh2); 