% Advanced script: Basis Pursuit using unconstrained optimization tool SesopTn
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

%  Michael  Zibulevsky 07.07.2008; 17.07.2008; 03.08.2008
%
% Copyright (c) 2008. All rights reserved. No warranty. Free for academic use


addpath('../');addpath('../util_mz');
addpath('../../minFunc')

%profile off
profile on
ColorPlots=0;  % 1 - use color plots; 0 - bw

INIT_DATA=1;         % 1 - init random data;  0 - use old data from the previous run
CONTINUE_OLD_RUN=0;  % 1- We can run again starting from  achieved solution with changed smoothing/concavity parameters
maxCG=[1 5 10 40];

options=sesoptn_optionset;     % Get default options structure (see comments in optionset_sesoptn.m) 
options.max_newton_iter=5;     % Max Newton iterations in subspace optimization
options.max_iter_CGinTN = 20;  % Conj.Grad steps in Truncated Newton (if  set to 0, TN  not used)
options.max_sesop_iter=1000/(options.max_iter_CGinTN+1);    % Max  SESOP iterations

options.FlagNonlinCG    = 0;  % 1 - Polak-Ribiere Nonlinear Conjugate Gradient method; 0 - any other method

options.nLastSteps=1;         % SESOP subspace will include: nLastSteps + [precond] gradient + [3 TN directions] + [2 Nemirovski directions])
options.PureTN=1;             %  1 - Pure Truncated Newton; 0 - Truncated Newton with  SESOP (recommended); 
options.FlagNemirovski_dir=0; % 1 - Add two directions of Nemirovski: x-x0 and  sum w_k*g_k

options.FlagRestrictedMatrix_TN=0; % 1 - restrict TN matrix; 0 - don't
options.PeriodRestoreAx=1;
options.PeriodRestoreAp=1;

options.dxTol=-1e-60;         % stopping criteria norm(change_in_x)
options.dfTol=-1e-60;         % abs(change_in_objective_func)
options.gTol=1e-6;            % norm(gradient)




options.sesop_figure_name=sprintf('SESOPtn  %d CG steps per TN iter',options.max_iter_CGinTN);

							   



n=400;   % measurement sample size
k=n;   % num of signal coefs to be recovered  = [dimensionality of the optimization problem]
%k=20*n;   % num of signal coefs to be recovered  = [dimensionality of the optimization problem]

sigma_noise=0.01;   % STD of white Gaussian noise added in signal generation
sparsity=  0.1;    % Sparsity of the generated coefs vector c00

par.weight_abs_penalty= 0; %0.01; % Weight of smooth abs penalty (mu)
par.eps_smooth_abs=1e-3;      % Smoothing parameter of asbsolute value approximation


if INIT_DATA    &  ~CONTINUE_OLD_RUN ,  
	%A=[hadamard_matrix(n) eye(n)];
	A=(1/sqrt(n))*randn(n,k);
	%A=spdiag(exp(-linspace(0,3,n)));
	%A= 0.01*A+diag(randn(n,1));  %To test diag. precond
end

par.multA= @(x,par)  multMatr(A,x);     % user function   y=Ax
par.multAt=@(x,par)  multMatrAdj(A,x);  % user function  y=A'*x

if INIT_DATA    &  ~CONTINUE_OLD_RUN ,  
	%c00=sign(sprandn(k,1,sparsity));    %  generate original vector of coefficients
	c00=1*sign(sprandn(k,1,sparsity));    %  generate original vector of coefficients
	x00=par.multA(c00,par);             % create clean measurement signal  
	y=x00(:)+sigma_noise*randn(n,1);    %  add noise 
	%c0=rand(size(c00));                 % Starting point for optimization
	c0=ones(size(c00));                 % Starting point for optimization
	par.x00=x00;
	par.y=y;

	% Compute diag_AtA  to be  used for preconditioning
    % diag_AtA=StochasticCalcDiagAtA(par.multAt,size(par.y),30,par); % For large matrices
	diag_AtA=diag(A'*A);                                       % For small matrices
end


par.func_u =@ls_fgh;                    % user function  f(u) = 0.5*||u - y||^2
par.func_x=@sum_fast_abs_smooth_new;   % user function  f(x) = mu*sum(abs_smoothed_eps(x))
%par.func_x=@sum_fast_abs_smooth;        % user function  f(x) = mu*sum(abs_smoothed_eps(x))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%      User preconditioning function (active, when options.precond=1):
%
%                   d_new=options.mult_precond(-g, x, par);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


options.mult_precond = @(g, x, Ax, InsideCGforTN, par) mult_diag_precond (g, x, Ax, InsideCGforTN, par,diag_AtA);





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%      Perform SESOP optimization
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




if CONTINUE_OLD_RUN,   c_init=c;
else                   c_init=c0;
end

par.multA= @(x,par)  multMatr(A,x);     % user function   y=Ax
par.multAt=@(x,par)  multMatrAdj(A,x);  % user function  y=A'*x

% options.max_iter_CGinTN=10;
% tic
% %c=sesoptn(c_init,par.func_u, par.func_x, par.multA, par.multAt,options,par);
% [c,report1]=sesoptn_svm(c_init,par.func_u, par.func_x, par.multA, par.multAt,options,par);
% toc
% 


if ColorPlots
	colorvec = [{'b'},{'r'},{'g'},{'y'},{'k'},{'c'},{'m'},{'.'},{'--'},{'r.'}];
else
	%colorvec = [{'k:'},{'k-.'},{'k--'},{'k.'},{'k'},{'k+'},{'k^'},{'ks'},{'k'}];
	%colorvec = [{'k'},{'k:'},{'k-.'},{'k--'},{'k.'},{'k+'},{'k^'},{'ko'},{'ks'}];
	colorvec = [{'k'},{'k:'},{'k-.'},{'ko'},{'k+'},{'k^'},{'ko'},{'ks'},{'k--'}];
	%colorvec = [{'k'},{'k.'},{'ko'},{'k*'},{'k^'},{'ks'},{'kd'},{'k>'},{'kv'}];
end



par.multA= @(x,par)  multMatr(A,x);     % user function   y=Ax
par.multAt=@(x,par)  multMatrAdj(A,x);  % user function  y=A'*x

kkk=length(maxCG);


j=0;
for i=1:kkk,
	options.max_iter_CGinTN=maxCG(i);
	options.max_sesop_iter=300/(options.max_iter_CGinTN+1);    
	[c,report(i)]=sesoptn(c_init,par.func_u, par.func_x, par.multA, par.multAt,options,par);
end	



fighandle1=figure('name', 'Function values with iterations'); 
fighandle2=figure('name', 'Function values with time'); 
legendstrings=[];

for i=1:kkk,
	j=j+1;
	ff=report(i).func_values;
	nncg=[0:length(ff)-1]*(maxCG(i) +1);
	figure(fighandle1);
	semilogy(nncg,ff,colorvec{i}); hold on
    %subplot(121);semilogy(nncg,ff,colorvec{i}); hold on
	%subplot(122);semilogy(nncg,report(i).gradnorms,colorvec{i}); hold on
	legendstrings{i}=sprintf('%d CG steps in TN',maxCG(i));

	figure(fighandle2);
	subplot(121);semilogy(report(i).times,ff,colorvec{i}); hold on
	subplot(122);semilogy(report(i).times,report(i).gradnorms,colorvec{i}); hold on

end
%figure(fighandle1);subplot(121);legend(legendstrings);%subplot(122);legend(legendstrings);
%figure(fighandle2);subplot(121);legend(legendstrings);%subplot(122);legend(legendstrings);
figure(fighandle1);legend(legendstrings);
xlabel('Global CG iteration count'); ylabel('Objective Function')

figure(fighandle2);legend(legendstrings);

figure;plot([c00+1,  c-1]); % add +-1 to the data in order to separate the plots
title('Basis Pursuit results:\newline Original (blue) and recovered (green) coefs');


%profile viewer;
return














%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Standard CG for quadratic problem min ||Ac-y||^2
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
par.flagXnew=1;
par.multA= @(x,par)  multMatr(A,x);     % user function   y=Ax
par.multAt=@(x,par)  multMatrAdj(A,x);  % user function  y=A'*x
optionsCG.max_iter_cg=2000;
optionsCG.period_show_progress=1;
multH=@(x,par) par.multAt(par.multA(x,par),par);
b=par.multAt(par.y,par);
c2=CGmz(multH,c_init,b,optionsCG,par); 
f_cg=fg(c2,par)

return




par.flagXnew=1;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% %
% %  Optimization using fminunc.m from Matlab Optimization toolbox 
% %
% %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% disp('Optimization using fminunc.m from Matlab Optimization toolbox')
% options1 =optimset('LargeScale','off','GradObj','on','Display','iter','MaxIter',200);
% c1 =fminunc( @fg, c_init, options1, par);






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%  Optimization using minFunc of Mark Schmidt 
%                              http://www.cs.ubc.ca/~schmidtm
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Optimization using minFunc of Mark Schmidt')
options1 =optimset('LargeScale','off','GradObj','on','Display','iter');
options1.Display='final';
options1.Method = 'lbfgs';
options1.TolFun=1e-16;
options1.TolX=1e-16;
options1.Corr=8;       % Number of previous steps used in L-BFGS update
options1.MaxFunEvals=4000;
options1.MaxIter=2000;

par.multA= @(x,par)  multMatr(A,x);     % user function   y=Ax
par.multAt=@(x,par)  multMatrAdj(A,x);  % user function  y=A'*x

tic;
[c1,f] =minFunc( @fg, c_init, options1,par);
toc
fprintf('f_opt=%g ',f)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%              Show results 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


figure;plot([c00+1,  c-1]); % add +-1 to the data in order to separate the plots
title('Basis Pursuit results:\newline Original (blue) and recovered (green) coefs');

profile viewer
