% UTIL_MZ: utilities
%
% Files
%   CGmz                    - Solve equation Hx = b with Conjugate Gradients, 
%   StochasticCalcDiagAtA   - Compute diag A'A: norms of rows of A'  as empirical variance
%   abs_smoothed_eps        - Convex:   |x| - eps* log(|x|/eps +1);  "Concave":  (1/h)*(log(1+c|x|) -  (1/p1)*log(1+p1*c|x|)); 
%   concave_lognorm         - Approximation of concave "norm"  phi(t)= (1/h)*(log(1+ct) -  (1/p)*log(1+pct))
%   concave_lognorm_preset  - Prepare parameters for the approximation of concave "norm" using logarithms
%   cubic_roots             - Solve a cubic equation   a*x^3 + b*x^2 + c*x + d = 0,     where a, b, c, and d are real.
%   diag_quadr_penalty      - f=0.5*par.quadrpenpar*sumsqr(x(:));
%   fast_abs_smooth_mz      - Smoothed approximation of absolute value function and its 1st and 2nd derivatives.  
%   fg                      - Value and gradient of f=par.f_u(Ax)+ par.f_x(x); For use by fminunc 
%   hadamard_matrix         - Create Hadamard matrix NxN
%   ls_fgh                  - Least-squares objective function:  f=||u-par.y||^2
%   muldm                   - Compute efficiently  B= diag(d)*A
%   mulmd                   - Compute efficiently  B= A*diag(d)
%   multMatr                - y=A*x
%   multMatr1               - y = Ax or A'x
%   multMatrAdj             - y=A'*x
%   mult_global_matrix      - Compute  y=GlobalMatrix*x
%   mult_global_matrix_adj  - Compute  y=GlobalMatrix' * x
%   spdiag                  - Diagonal matrix in sparse format
%   sum_abs_smooth          - Smoothed L1-norm penalty:  mu*sum(abs_smooth(x));  
%   sum_fast_abs_smooth_new - Penalty for  vector x:  mu*sum(|x|), 
