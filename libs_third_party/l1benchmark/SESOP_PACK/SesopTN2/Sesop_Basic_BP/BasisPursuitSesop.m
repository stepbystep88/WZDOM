% Simple script: Basis Pursuit using unconstrained optimization tool SesopTn
%
% Find sparse solution of underdetermined linear system of equations Ac=y using L1-norm penalty
%
%  A -  "Wide" system matrix (dictionary)
%  y = Ac + noise  - measurement vector
%  c -  vector coefficients to be recovered using optimization
%
%  min_c  ||Ac-y||^2 +  mu*sum(smooth_abs(c))
%
% where smooth_abs(c) - smoothed approximation asbsolute value function 
%

%  Michael  Zibulevsky 07.07.2008; 17.07.2008; 03.08.2008; 25.02.2010
% Copyright (c) 2008. All rights reserved. No warranty. Free for academic use

%profile off
%profile on
addpath('../');addpath('../util_mz');

INIT_DATA=1;         % 1 - init random data;  0 - use old data from the previous run
CONTINUE_OLD_RUN=0;  % 1- We can run again starting from  achieved solution with changed smoothing/concavity parameters

options=sesoptn_optionset;   % Get default options structure (see comments in optionset_sesoptn.m) 
options.max_sesop_iter=200;  % Max  SESOP iterations
options.max_newton_iter=1;   % Max Newton iterations in subspace optimization (one can play with this)
options.max_iter_CGinTN = 0; %10;  % Conj.Grad steps in Truncated Newton (if  set to 0, TN  not used)

options.sesop_figure_name=sprintf('SESOPtn  %d CG steps per TN iter',options.max_iter_CGinTN);

n=128;   % measurement sample size
k=2*n;   % num of signal coefs to be recovered  = [dimension of the optimization problem]

sigma_noise=0.01;   % STD of white Gaussian noise added in signal generation
sparsity=  0.03;    % Sparsity of the generated coefs vector c00

par.weight_abs_penalty= 1e-2; % Weight of smooth abs penalty (mu)
par.eps_smooth_abs=1e-2;      % Smoothing parameter of asbsolute value approximation


if INIT_DATA    &  ~CONTINUE_OLD_RUN ,  
	A=[hadamard_matrix(n) eye(n)];
	%A=(1/sqrt(n))*randn(n,k);
	%A= 0.01*A+diag(randn(n,1));  %To test diag. precond
end

par.multA= @(x,par)  multMatr(A,x);     % user function   y=Ax
par.multAt=@(x,par)  multMatrAdj(A,x);  % user function  y=A'*x

if INIT_DATA    &  ~CONTINUE_OLD_RUN ,  
	c00=sign(sprandn(k,1,sparsity));    %  generate original vector of coefficients
	x00=par.multA(c00,par);             % create clean measurement signal  
	y=x00(:)+sigma_noise*randn(n,1);    %  add noise 
	c0=rand(size(c00));                 % Starting point for optimization
	par.x00=x00;
	par.y=y;

	% Compute diag_AtA  to be  used for preconditioning
   % par.diag_AtA=StochasticCalcDiagAtA(par.multAt,size(par.y),30,par); % For large matrices
	par.diag_AtA=diag(A'*A);                                       % For small matrices
end


par.func_u =@ls_fgh;                   % user function  f(u) = 0.5*||u - y||^2
par.func_x=@sum_fast_abs_smooth_new;   % user function  f(x) = mu*sum(abs_smoothed_eps(x))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%      User preconditioning function (active, when options.precond=1):
%
%                   d_new=options.mult_precond(-g, x, par);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% options.precond=1;  % uncomment these two lines when using preconditioning
% options.mult_precond = @(g, x, Ax, InsideCGforTN, par) mult_diag_precond (g, x, Ax, InsideCGforTN, par);





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%         Perform SESOP optimization
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




if CONTINUE_OLD_RUN,   c_init=c;
else                   c_init=c0;
end


c=sesoptn(c_init,par.func_u, par.func_x, par.multA, par.multAt,options,par);





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%              Show results 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


subplot(224); plot([c00+1,  c-1]); % add +-1 to the data in order to separate the plots
title('Basis Pursuit results:\newline Original (blue) and recovered (green) coefs');

%profile viewer
