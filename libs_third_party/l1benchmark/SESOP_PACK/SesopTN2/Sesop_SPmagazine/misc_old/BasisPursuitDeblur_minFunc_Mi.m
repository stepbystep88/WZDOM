% Basis Pursuit using unconstrained optimization tool SesopTn
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

%  Michael  Zibulevsky 07.07.2008; 17.07.2008; 03.08.2008; 10.06.2009
%
% Copyright (c) 2008. All rights reserved. No warranty. Free for academic use


%addpath('..');
addpath('../');addpath('../util_mz');
addpath('../../minFunc')
addpath('../../Rice_Wavelab')

global GlobalNiter GlobalGradNorms GlobalFuncValues GlobalNNiter GlobalXsnr; 
GlobalNNiter=[];GlobalNiter=0;GlobalGradNorms=[]; GlobalFuncValues=[]; GlobalXsnr=[];

%profile off
%profile on

%
INIT_DATA=0;         % 1 - init random data;  0 - use old data from the previous run
CONTINUE_OLD_RUN=0;  % 1- We can run again starting from  achieved solution with changed smoothing/concavity parameters

FlagMinFunc = 0;     % Optimization using L-BFGS, minFunc of Mark Schmidt
FlagSESOP   = 1;     % Optimization using SESOP

options=sesoptn_optionset;  % Get default options structure (see comments in optionset_sesoptn.m) 
options.max_sesop_iter  = 60;  % Max  SESOP iterations
options.max_newton_iter = 7;    % Max Newton iterations in subspace optimization
options.max_iter_CGinTN = 0;   % Conj.Grad steps in Truncated Newton (if  set to 0, TN  not used)

options.precond = 1;         % 1 - use user-defined  preconditioner, 0 - don't
FlagPCD         = 1;

options.period_show_progress=1;    % Periodicity of showing sesop and user plots


options.report_func=@deblur_report_func;   % User function reference to display/collect iteration progress;
options.sesop_figure_name=sprintf('SESOPtn  %d CGstepsPerTN; precond=%d PCD=%d',options.max_iter_CGinTN,options.precond, FlagPCD);

par.report_func=options.report_func;
par.period_show_progress=options.period_show_progress;


%k=4*n;   % num of signal coefs to be recovered  = [dimensionality of the optimization problem]
n=128;    % signal/image size
[ii,jj]=meshgrid(-7:1:7); 
h=1./(ii.^2+jj.^2+1); 
h=h/sum(h(:)); 
% h=gausswin(7);      % prepare for blur kernel
sigma_noise=4;   % STD of white Gaussian noise added in signal generation
%sparsity=  0.03;   % Sparsity of the generated coefs vector c00

par.weight_abs_penalty= 0.2; % Weight of smooth abs penalty (mu)
par.eps_smooth_abs=1e-2;      % Smoothing parameter of asbsolute value approximation



par.wavelet_cqf = daubcqf(4,'min');  % Wavelet definitions
par.wavelet_levels = 3;              %
par.SynthesisFunc = @TI2Dsynthesismz;% 
par.AnalysisFunc  = @TI2Danalysismz; % 



par.blur_kernel=h ; % conv2(h,h');          %figure;mesh(h2);title('Blur kernel');
par.projector    =@blurring;
par.projector_adj=@blurring_adj;

% par.projector    = @(x,par) x;  % Identity projector
% par.projector_adj= @(x,par) x;

par.multA=@projector_synthesis;
par.multAt=@analysis_projector_adj;



%
if INIT_DATA    &&  ~CONTINUE_OLD_RUN ,  
   X00=256*phantom(n);          % Create clean image
   par.imagesize=size(X00);
   par.proj_size=par.imagesize;
   [coefs,coefs_size]=feval(par.AnalysisFunc, X00(:),par);
   %[coefs,coefs_size]=par.AnalysisFunc(X00(:),par);      % this line does not work in matlab2004
   par.coefs_size=coefs_size;
   
   %c00=sign(sprandn(k,1,sparsity));    %  generate original vector of coefficients
	%x00=par.multA(c00,par);             % create clean measurement signal  
	
	par.x00=X00(:);
   par.y00=par.projector(X00(:),par);
   par.y=par.y00+sigma_noise*randn(size(par.y00));    %  add noise 
   
	c0=par.multAt(par.y,par);             % Starting point for optimization

	fprintf('\n Compute diag_AtA  to be  used for preconditioning...')
	diag_AtA=StochasticCalcDiagAtA(par.multAt,size(par.y),300,par); % For large matrices
	%diag_AtA=diag(A'*A);                                       % For small matrices
   fprintf('Done\n');
   figure;plot(diag_AtA);title('diag_AtA')
   
%    fprintf('\n Compute max singular value of A used in SSF...');
%    sigma_max=max_singular_value(par.multA,par.multAt,ones(size(c0)),50,par)
%    fprintf('Done\n');
end


par.func_u =@ls_fgh;                    % user function  f(u) = 0.5*||u - y||^2
par.func_x=@sum_abs_smooth;   % user function  f(x) = mu*sum(abs_smoothed_eps(x))
%par.func_x=@sum_fast_abs_smooth_new;   % user function  f(x) = mu*sum(abs_smoothed_eps(x))


%
%
%   User preconditioning function: d_new=mult_precond (-g, x, par)
%   (active, when options.precond=1)
%


if FlagPCD   % PCD direction 

	options.mult_precond =@(g, x, Ax, InsideCGforTN, par) mult_precond_pcd(g, x, Ax, InsideCGforTN, par,diag_AtA);

else         % Diagonal preconditioning
	
	options.mult_precond = @(g, x, Ax, InsideCGforTN, par) mult_diag_precond(g, x, Ax, InsideCGforTN, par,diag_AtA);

end


%options.mult_precond = @(g, x, Ax, InsideCGforTN, par) mult_diag_precond (g, x, Ax, InsideCGforTN, par,diag_AtA);





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

if FlagSESOP,
   par.BPblur_report_figure_handle=figure('Position',[10 40 400 600],'name','SESOP BP_blur_report');
   tic
   c=sesoptn(c_init,par.func_u, par.func_x, par.multA, par.multAt,options,par);
   toc
end
%profile viewer;
%return

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% % Standard CG for quadratic problem min ||Ac-y||^2
% %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% par.multA= @(x,par)  multMatr(A,x);     % user function   y=Ax
% par.multAt=@(x,par)  multMatrAdj(A,x);  % user function  y=A'*x
% optionsCG.max_iter_cg=2000;
% optionsCG.period_show_progress=10;
% multH=@(x,par) par.multAt(par.multA(x,par),par);
% b=par.multAt(par.y,par);
% c2=CGmz(multH,c_init,b,optionsCG,par); 
% f_cg=fg(c2,par)
% 
% return


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

options1 =optimset('GradObj','on','Display','iter');
%options1.Display='final';
options1.Method = 'lbfgs';
options1.TolFun=1e-16;
options1.TolX=1e-16;
options1.Corr=8;       % Number of previous steps used in L-BFGS update
options1.MaxFunEvals=4000;
options1.MaxIter=options.max_sesop_iter; %2000;

%par.multA= @(x,par)  multMatr(A,x);     % user function   y=Ax
%par.multAt=@(x,par)  multMatrAdj(A,x);  % user function  y=A'*x

if FlagMinFunc
   disp('Optimization using minFunc of Mark Schmidt')
   par.BPblur_report_figure_handle=figure('Position',[10 40 400 600],'name','L-BFGS BP_blur_report');
   tic;
   [c1,f] =minFunc( @fg, c_init, options1,par);
   toc
   fprintf('f_opt=%g ',f)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%              Show results 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%subplot(224); plot([c00+1,  c-1]); % add +-1 to the data in order to separate the plots
%title('Basis Pursuit results:\newline Original (blue) and recovered (green) coefs');

%profile viewer
