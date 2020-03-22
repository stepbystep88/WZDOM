% "Concave norm" Basis Pursuit using unconstrained optimization tool SesopTn and PCD (parallel coordinate descent)
%
% Find sparse solution of underdetermined linear system of equations Ac=y
% using approximation of Lp-norm penalty,  0<p<2
%
%  A -  "Wide" system matrix (dictionary)
%  y = Ac + noise  - measurement vector
%  c -  vector coefficients to be recovered using optimization
%
%  min_c  ||Ac-y||^2 +  mu*sum(smooth_abs(c))
%
% where smooth_abs(c) - smoothed approximation of |c|^p,  0<p<2
%

%  Michael  Zibulevsky 07.07.2008; 17.07.2008; 04.08.2008; 04.02.2009
%
% Copyright (c) 2008 - 2009. All rights reserved. No warranty. Free for academic use

%profile off
%profile on
addpath('../');addpath('../util_mz');

global Global_conc_lognorm;

INIT_DATA=1;         % 1 - init random data;  0 - use old data from the previous run
CONTINUE_OLD_RUN=0;  % 1- We can run again starting from  achieved solution with changed smoothing/concavity parameters

options=sesoptn_optionset;  % Get default options structure (see comments in optionset_sesoptn.m) 
options.max_iter_CGinTN=0;  % Conj.Grad steps in Truncated Newton (if  set to 0, TN  not used)

options.precond=1;      % 1 - use user defined  preconditioning,  0 - don't use preconditioning

FlagPCD=1;              % when options.precond=1:
                        %     1 -  PCD (parallel coord. descent)
                        %     0 - diagonal precond.
						
Global_conc_lognorm.on = 0;      % 1 - use Lp_norm smooth approximation, p<1, using two log's
%                                  0 - use convex smooth |c| approximation

% It is recommended to perform initial optimization with convex |c| approximation: Global_conc_lognorm.on = 0;
% then continue from the achieved point with concave one (uncomment the following line):
%
%Global_conc_lognorm.on = 1; CONTINUE_OLD_RUN=1;


options.sesop_figure_name=sprintf('SESOPtn  %d CG steps per TN iter; FlagPCD=%d',options.max_iter_CGinTN, FlagPCD);

n=128;   % measurement sample size
k=2*n;   % num of signal coefs to be recovered  = [dimensionality of the optimization problem]

sigma_noise=0.01;   % STD of white Gaussian noise added in signal generation
sparsity=  0.03;    % Sparsity of the generated coefs vector c00

par.weight_abs_penalty= 1e-2; % Weight of smooth abs penalty (mu)
par.eps_smooth_abs=1e-4;      % Smoothing parameter of asbsolute value approximation





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%           Using "concave" norm approximation in the second term of 
%
%                   min_c  ||Ac-y||^2 +  mu*sum(smooth_abs(c))
%
%    It is recommended to minimize convex problem first, and then minimize
%    concave function, setting       Global_conc_lognorm.on = 1
%    using achieved starting point:  CONTINUE_OLD_RUN    = 1
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% we use phi(t)= (1/h)*(log(1+ct) -  log(1+pct)/p)
% where  h=  (log(1+c) - log(1+p*c)/p)
	
c_conc=3    % argument scaling
b_conc=100  %1./par.eps_smooth_abs; % second derivative at the origin
                                   % (may be automatically increased by concave_lognorm_preset to match c)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%        SESOP stopping criteria: if any of them is satisfied, SESOP terminates
options.dxTol = 1e-6;    % norm(change_in_x)
options.dfTol = 1e-12;   % abs(change_in_objective_func)
options.gTol  = 1e-6;    % norm(gradient)



																		  
if INIT_DATA    &  ~CONTINUE_OLD_RUN ,  
	A=[hadamard_matrix(n) eye(n)];
	%A=(1/sqrt(n))*randn(n,k);   
	%A= 0.01*A+diag(randn(n,1));  %To test diag. precond
end

par.multA= @(x,par)  multMatr(A,x);        % user function   y=Ax
par.multAt=@(x,par)  multMatrAdj(A,x);     % user function  y=A'*x

if INIT_DATA    &  ~CONTINUE_OLD_RUN ,  
	c00=sign(sprandn(k,1,sparsity));    %  generate original vector of coefficients
	x00=par.multA(c00,par);             % create clean measurement signal  
	y=x00(:)+sigma_noise*randn(n,1);    %  add noise 
	c0=rand(size(c00));                 % Starting point for optimization
	par.x00=x00;
	par.y=y;

	% Compute diag_AtA  to be  used for preconditioning
    % diag_AtA=StochasticCalcDiagAtA(par.multAt,size(par.y),20,par); % For large matrices
	diag_AtA=diag(A'*A);   % For small matrices
end


par.func_u =@ls_fgh;         % user function  f(u) = 0.5*||u - y||^2
par.func_x=@sum_abs_smooth;  % user function  f(x) = mu*sum(abs_smoothed_eps(x))


if Global_conc_lognorm.on    % Lp_norm smooth approximation, p<1, using two log's
	%
	% h=  (log(1+c) - log(1+p*c)/p)
	%
	[p_tmp, h_tmp, b_tmp] = concave_lognorm_preset(b_conc,  c_conc)
	% We will use the calculated parameters to build the function concave_lognorm:
	% phi(t)= (1/h)*(log(1+ct) -  log(1+pct)/p)
end






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%      User preconditioning function: d_new=mult_precond (-g, x, par);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




if FlagPCD   % PCD "preconditioning"

	options.mult_precond =@(g, x, Ax, InsideCGforTN, par) mult_precond_pcd(g, x, Ax, InsideCGforTN, par,diag_AtA);

else         % Diagonal preconditioning
	
	options.mult_precond = @(g, x, Ax, InsideCGforTN, par) mult_diag_precond (g, x, Ax, InsideCGforTN, par,diag_AtA);

end







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%         Perform SESOP optimization
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




if CONTINUE_OLD_RUN,     c_init=c;
else                     c_init=c0;
end

c=sesoptn(c_init,par.func_u, par.func_x, par.multA, par.multAt,options,par);

% Try without preconditioning:
 % options.precond=0;
 % c=sesoptn(c0,par.func_u, par.func_x, par.multA, par.multAt,options,par);

 
 
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%              Show results 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


subplot(224); plot([c00+1,  c-1]); % add +-1 to the data in order to separate the plots
title('Basis Pursuit results:\newline Original (blue) and recovered (green) coefs');

%profile viewer
