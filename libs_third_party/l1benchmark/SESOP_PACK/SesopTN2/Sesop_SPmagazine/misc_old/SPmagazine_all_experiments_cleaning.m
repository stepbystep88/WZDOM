% Deblur, Tomography, Compr. Sensing, "bad" matrices
%
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

%  Michael  Zibulevsky 07.07.2008; 17.07.2008; 03.08.2008;
%  10.06.2009; 04.08.2009; 23.08.2009, 02.2010
%
% Copyright (c) 2008. All rights reserved. No warranty. Free for academic use

%addpath('..');
addpath('../');addpath('../util_mz'); addpath('../tomogr_regul');
addpath('../../minFunc')
addpath('../../Rice_Wavelab')
addpath('../../l1_ls_matlab')

global GlobalNiter GlobalGradNorms GlobalFuncValues GlobalTimes GlobalNNiter GlobalNNiterFG GlobalSNRtime GlobalXsnr GlobalY_RNR; 
global          Global_conc_lognorm; 
Global_conc_lognorm.simple=1;


flagSESOP_TN_10=0;
flagPCD_SESOP_TN_10=0;
flagSESOP_7=0;
flagSESOP_1=0;
flagSESOP_0=0;
flagCG=0;
flagPCD=0;
flagPCD_CG=0;
flagPCD_SESOP_7=0;
flagSSF=0;
flagSSF_lsrch=0;
flagSSF_CG=0;
flagSSF_SESOP_7=0;
flagLBFGS=0;
flagL1_LS_IntPoint=0;
flagNesMz_SSF=0;
flagFISTA_SSF=0;
flagPRE_CG=0;
flagPRE_SESOP=0;
flagSesopLBFGS=0;
flagSSF_LBFGS=0;
flagPCD_LBFGS=0;


% flagPCD=1;
 flagPCD_SESOP_7=1;
% flagSSF=1;
% flagSSF_SESOP_7=1;
% flagFISTA_SSF=1;



maxNumCompThreads(1); % To avoid multiple CPU use (for correct CPU time measuring)
FlagDenoise=0;FlagSparseImage=0; image_sparsity=1;resize_factor=1;


options=sesoptn_optionset;  % Get default options structure (see comments in optionset_sesoptn.m) 
par.ProblSubName='';

INIT_DATA=1;                % 1 - init random data;  0 - use old data from the previous run

n0=64;  m0=n0              % Image size: m0-by-n0 


par.ImageName='Phantom';
%par.ImageName='Sparse'; FlagSparseImage=1; image_sparsity=0.05; 
%par.ImageName='Lena'; 
%par.ImageName='Peppers'; 



par.ProblName='Deblur'; FlagDenoise=0;   % when FlagDenoise=1 - use delta-function blur (just denoising problem)
%par.ProblName='ComprSens';
%par.ProblName='Tomogr';
%par.ProblName='Matrix';par.ImageName='Sparse'; FlagSparseImage=1; image_sparsity=0.05; 


%matrdir='C:\Users\user\Videos\timingdata\timingdata\finaldata\';
%par.ProblSubName='badgaussmatrix'; A=load([matrdir 'badgaussmatrix']);A=A.mat; 
%par.ProblSubName='tonytaperwavreal';A=load([matrdir 'tonytaperwavreal']); A=A.mat; A = (1/2.0472e-004)*A;
%par.ProblSubName= 'geogaussmatrix'; A=load([matrdir 'geogaussmatrix']); A=A.mat; A = (1/2.0472e-004)*A;
%par.ProblSubName= 'gaussmatrix'; A=load([matrdir 'gaussmatrix']); A=A.mat;




par.weight_abs_penalty                 = 1e-6; % Weight of smooth abs penalty (mu)
%par.weight_abs_penalty                 = 1e-3; % Weight of smooth abs penalty (mu)

if FlagDenoise, par.weight_abs_penalty = 0.4e-1; end 

%options.max_sesop_iter  = 50;    % Max  SESOP iterations
%options.max_sesop_iter  = 1300;    % Max  SESOP iterations
options.period_show_progress=5;  % Periodicity of showing sesop and user plots

%options.max_sesop_iter  = 200;    % Max  SESOP iterations
%options.period_show_progress=20;  % Periodicity of showing sesop and user plots

if  strcmp(par.ProblName,'ComprSens')
   options.max_sesop_iter  = 1000;    % Max  SESOP iterations
   options.period_show_progress=100;  % Periodicity of showing sesop and user plots
   par.weight_abs_penalty                 = 1e-3; 
end

% options.max_sesop_iter  = 1500;    % Max  SESOP iterations
% options.period_show_progress=100;  % Periodicity of showing sesop and user plots


FlagConcaveNorm=0;                  % 1 - use concave norn; 0 - convex one


options.ContinueOldRun=0;           % 1- We can run again starting from  achieved solution with changed smoothing/concavity parameters
if INIT_DATA, options.ContinueOldRun=0;end

 
options.max_newton_iter = 7; % Max Newton iterations in subspace optimization
options.nLastSteps=1;        % SESOP subspace will include: nLastSteps + [[precond] gradient] +

if  strcmp(par.ProblName,'Matrix')
   par.sigma_noise                = 1e-4;   % STD of white Gaussian noise in measurements
else
   par.sigma_noise                = 1e-2;   % STD of white Gaussian noise in measurements
end
if FlagDenoise,par.sigma_noise = 6e-2; end



par.eps_smooth_abs=1e-3;      % Smoothing parameter of asbsolute value approximation
%par.eps_smooth_abs=1e-8;      % Smoothing parameter of asbsolute value approximation



Global_conc_lognorm.on = 0;      % 1 - use Lp_norm smooth approximation, p<1, using two log's
%                                  0 - use convex smooth |c| approximation

% It is recommended to perform initial optimization with convex |c| approximation: Global_conc_lognorm.on = 0;
% then continue from the achieved point with concave one (uncomment the following line):
%
%Global_conc_lognorm.on = 1; CONTINUE_OLD_RUN=1;




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
	
c_conc=3;    % argument scaling
b_conc=100;  %1./par.eps_smooth_abs; % second derivative at the origin
                                   % (may be automatically increased by concave_lognorm_preset to match c)

                                   
% if Global_conc_lognorm.on    % Lp_norm smooth approximation, p<1, using two log's
	%
	% h=  (log(1+c) - log(1+p*c)/p)
	%
   [p_tmp, h_tmp, b_tmp] = concave_lognorm_preset(b_conc,  c_conc);
	% We will use the calculated parameters to build the function concave_lognorm:
	% phi(t)= (1/h)*(log(1+ct) -  log(1+pct)/p)
% end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%











options.PeriodRestoreAx=8;    % period of SESOP restorations to avoid error accumulation
options.PeriodRestoreAp=1;

options.ShowSesopPlots=1;           % 1 - plot iteration progress, 0 - don't


options.dxTol = 1e-16;    % norm of change_in_x
options.dfTol = 1e-16;   % abs change_in_objective_func
options.gTol  = 1e-16;    % norm of the gradient


%options.max_iter_CGinTN = 0;   % Conj.Grad steps in Truncated Newton (if  set to 0, TN  not used)
%options.max_sesop_iter=ceil(options.max_sesop_iter/(options.max_iter_CGinTN+1));

options.report_func=@deblur_report_func;   % User function reference to display/collect iteration progress;
options.sesop_figure_name=sprintf('SESOPtn'); %d CGstepsPerTN; precond=%d PCD=%d SSF=%d', options.max_iter_CGinTN,options.precond, par.FlagPCD,par.FlagSSF);

par.report_func=options.report_func;
par.period_show_progress=options.period_show_progress;


%k=4*n;   % num of signal coefs to be recovered  = [dimensionality of the optimization problem]
%sparsity=  0.03;   % Sparsity of the generated coefs vector c00

if FlagSparseImage,
   par.SynthesisFunc = @(c,par) c;
   par.AnalysisFunc  = @eye_analysis;
else
   par.wavelet_cqf = daubcqf(4,'min');  % Wavelet definitions
   par.wavelet_levels = 4;              %
   par.SynthesisFunc = @TI2Dsynthesismz;%
   par.AnalysisFunc  = @TI2Danalysismz; %
end


%options.TESTGRAD=1;  % 1 - Test gradient numerically, plot testing results (without run optimization); 0 - don't test


if INIT_DATA 
   
   if strcmp(par.ProblName,'Matrix')
      %A=load([matrdir MatrName]); A=A.mat; 
      m0=64;n0=128;
      %X00=10*(sprandn(m0,n0,image_sparsity));save X00_sparse64_128 X00;
      par.projector     = @(x,par) A*x(:);
      par.projector_adj = @(y,par) (y(:)'*A)';
      par.proj_size=[8,231];
   end
   
   if FlagSparseImage,
      %X00=10*(sprandn(m0,n0,image_sparsity));%save X00 X00;
      if m0==64 && n0==128, load X00_sparse64_128; 
      else                  load X00_sparse64_64; 
      end
      X00=full(X00);
   else 
      if strcmp(par.ImageName,'Phantom'), X00=phantom(n0);end  % Create clean image
      if strcmp(par.ImageName,'Lena'),
         X00=imread('..\..\..\images\lena.png');
      end
      if strcmp(par.ImageName,'Peppers')
         X00=imread('..\..\..\images\peppers.png');
      end
      X00=double(X00);X00=X00/max(X00(:));%figure;imagesc(X00);colorbar
   end
   [m,n]=size(X00);  resize_factor=1/ceil(m/m0);
   if resize_factor~=1, X00=imresize(X00,resize_factor);  [m,n]=size(X00);end
   par.imagesize=[m,n];
   
   if strcmp(par.ProblName,'Deblur')
      %h=gausswin(7); h2=conv2(h,h');     % prepare for blur kernel
      [ii,jj]=meshgrid(-7:1:7); h2=1./(ii.^2+jj.^2+1);
      if FlagDenoise, h2=1; end % Delta-function blur (pure denoising problem)
      h2=h2/sum(h2(:));
      par.blur_kernel=h2; %figure;mesh(h2);title('Blur kernel');
      par.proj_size=par.imagesize;

   elseif strcmp(par.ProblName,'ComprSens')

      %nx=m*n;           % num of elements in image X
      %np=ceil(nx/4);    % num of projection samples
      %par.proj_size=[np,1];
      %P=1/n*randn(np,nx);  % Random projection matrix
      
      mlow=ceil(m/6); % Preserved low frequensy FFT components 
      nlow=ceil(n/6);
      %par.ComprSens.mlow=mlow; 
      %par.ComprSens.nlow=nlow;
      if FlagSparseImage,
         tmp=sprand(m,n,0.1);
      else
         tmp=sprand(m,n,0.05);
         tmp(1:mlow,1:nlow)=1;
         tmp(1:mlow,end-nlow:end)=1;
      end
      par.ComprSens.ind=find(tmp);
      par.proj_size=[2*numel(par.ComprSens.ind),1];

%      par.proj_size=[ceil(n),n];
%      ind=find(sprand(n,n,0.5));
%      P=1/sqrt(n)*randn(par.proj_size);  % Random projection matrices
%      Q=1/sqrt(n)*randn(par.imagesize)';
   end
end


if strcmp(par.ProblName,'Deblur')
   par.projector    =@blurring;
   par.projector_adj=@blurring_adj;

elseif strcmp(par.ProblName,'ComprSens')

   %par.projector     =@(x,par) P*x;
   %par.projector_adj =@(x,par) (x'*P)';

   par.projector     =@ComprSensProj;
   par.projector_adj =@ComprSensProj_adj;
   %par.projector     =@(x,par) reshape(P*reshape(x,par.imagesize)*Q, prod(par.proj_size),1)  ;
   %par.projector_adj =@(x,par) reshape(P'*reshape(x,par.proj_size)*Q', prod(par.imagesize),1)  ;

elseif strcmp(par.ProblName,'Tomogr'),
   par.flag_fastradon=1;
   par.k_small_bins=4;
   par.res=m;
   Nangles=2*par.res;
   par.angles = linspace(-45,135,Nangles); % Our fast Radon works only for angles from -45 to 135 degrees
   y=(1/n)*myradon(X00,par);
   %figure; subplot(121);imagesc(X00);subplot(122);imagesc(y);colorbar

   [Nbins,Nang]=size(y);
   par.Nbins=Nbins;
   par.Nang=Nang;
   par.proj_size=size(y);
   par.projector     = @(x,par) (1/n)*reshape(   myradon(reshape(x,par.imagesize),par), prod(par.proj_size),1);
   par.projector_adj = @(y,par) (1/n)*reshape(myadjradon(reshape(y,par.proj_size),par), prod(par.imagesize),1);
end

% par.projector    = @(x,par) x;  % Identity projector
% par.projector_adj= @(x,par) x;

par.multA=@projector_synthesis;
par.multAt=@analysis_projector_adj;



if INIT_DATA    &&  ~options.ContinueOldRun ,  
   
   [coefs,coefs_size]=feval(par.AnalysisFunc, X00(:),par);
   %[coefs,coefs_size]=par.AnalysisFunc(X00(:),par);      % this line does not work in matlab2004
   par.coefs_size=coefs_size;
   
   %c00=sign(sprandn(k,1,sparsity));    %  generate original vector of coefficients
	%x00=par.multA(c00,par);             % create clean measurement signal  
	
	par.x00=X00(:);
   par.y00=par.projector(X00(:),par);
   par.y=par.y00+par.sigma_noise*randn(size(par.y00));    %  add noise 
   
	fprintf('\n Compute diag_AtA  to be  used for preconditioning...')
	par.diag_AtA=StochasticCalcDiagAtA(par.multAt,size(par.y),100,par); % For large matrices
	%par.diag_AtA=diag(A'*A);                                       % For small matrices
   fprintf('Done\n');
   figure;plot(par.diag_AtA);title('diag_AtA')
   
   %c0=0*par.multAt(par.y,par);   
   c0=par.multAt(par.y,par);   
   
   fprintf('\n Compute max singular value of A used in SSF...');
   max_sing_value_A = max_singular_value(par.multA,par.multAt,randn(size(c0)),100,par)
   par.max_sing_value_A=max_sing_value_A;
   fprintf('Done\n');

	c0=(1/max_sing_value_A^2)*c0;    % Starting point for optimization

end

c0=par.multAt(par.y,par);
c0=(1/max_sing_value_A^2)*c0;    % Starting point for optimization


par.func_u =@ls_fgh;                    % user function  f(u) = 0.5*||u - y||^2
par.func_x=@sum_abs_smooth;   % user function  f(x) = mu*sum(abs_smoothed_eps(x))
%par.func_x=@sum_fast_abs_smooth_new;   % user function  f(x) = mu*sum(abs_smoothed_eps(x))


FlagMinFunc = 0;     % Optimization using L-BFGS, minFunc of Mark Schmidt
FlagSESOP   = 1;     % Optimization using SESOP
FlagBoyd    = 0;     % Optimization using Boyd's l1_ls interior point solver

options.precond = 1;         % 1 - use user-defined  preconditioner, 0 - don't
par.FlagPCD     = 0;
par.FlagSSF     = 0;
par.FlagPCDSSF_inside_TN=0;
options.FlagNonlinCG = 0;    % 1 - Polak-Ribiere Nonlinear Conjugate Gradient method; 0 - any other method
options.max_iter_CGinTN = 0; 

%f_figure_handle=figure('Position',[400 40 600 400],'name',['Objective function' par.ProblName] );


i_run=0;clear reports; 



if FlagConcaveNorm
   Global_conc_lognorm.on = 0;      % 1 - use Lp_norm smooth approximation, p<1, using two log's

   disp('Find  "warm start"  by convex optimization')
   method='SSF-CG-convex',
   options.precond = 1; FlagSESOP   = 1; par.FlagPCD=1;par.FlagSSF = 0;
   options.FlagNonlinCG = 1;
   options.nLastSteps=1;
   OptimizeBP; % Perform Optimization
   options.FlagNonlinCG = 0;

   disp('Concave penalty optimization')
   Global_conc_lognorm.on = 1,   par.eps_smooth_abs=1e0;   % 1 - use Lp_norm smooth approximation, p<1
   c0=c;  % Use "warm start" given by convex optimization
end


if flagSESOP_TN_10
i_run=i_run+1;method='SESOP-TN-10', if FlagConcaveNorm, method=[method '-concave'],end
FlagSESOP=1; par.FlagPCD=0;par.FlagSSF=0;FlagMinFunc = 0;
options.max_iter_CGinTN = 10; 
period_show_progress_old=options.period_show_progress;
%options.period_show_progress=ceil(options.period_show_progress/options.max_iter_CGinTN);
options.max_sesop_iter_old=options.max_sesop_iter;
options.max_sesop_iter=2*ceil(options.max_sesop_iter/options.max_iter_CGinTN);
options.precond = 0;         % 1 - use user-defined  preconditioner, 0 - don't
options.nLastSteps=5; 
OptimizeBP;  % Perform Optimization
Report.method=method;   reports(i_run)=Report;
options.max_iter_CGinTN = 0; 
options.period_show_progress=period_show_progress_old;
options.max_sesop_iter=options.max_sesop_iter_old;
end

if flagPCD_SESOP_TN_10
i_run=i_run+1;method='PCD-SESOP-TN-10', if FlagConcaveNorm, method=[method '-concave'],end 
FlagSESOP=1; par.FlagPCD=1;par.FlagSSF=0;FlagMinFunc = 0;
options.max_iter_CGinTN = 10; 
period_show_progress_old=options.period_show_progress;
%options.period_show_progress=ceil(options.period_show_progress/options.max_iter_CGinTN);
options.max_sesop_iter_old=options.max_sesop_iter;
options.max_sesop_iter=2*ceil(options.max_sesop_iter/options.max_iter_CGinTN);
options.nLastSteps=10; 
par.FlagPCDSSF_inside_TN=0;
OptimizeBP;  % Perform Optimization
Report.method=method;   reports(i_run)=Report;
options.max_iter_CGinTN = 0; 
options.period_show_progress=period_show_progress_old;
options.max_sesop_iter=options.max_sesop_iter_old;
end



if flagSESOP_7
i_run=i_run+1;method='SESOP-7', if FlagConcaveNorm, method=[method '-concave'],end
options.precond = 0;         % 1 - use user-defined  preconditioner, 0 - don't
FlagSESOP=1; par.FlagPCD=0;par.FlagSSF=0;FlagMinFunc = 0;
options.nLastSteps=7; 
OptimizeBP;  % Perform Optimization
Report.method=method;   reports(i_run)=Report;
end


if flagSESOP_0
i_run=i_run+1;method='SESOP-0', if FlagConcaveNorm, method=[method '-concave'],end
options.precond = 0;         % 1 - use user-defined  preconditioner, 0 - don't
FlagSESOP=1; par.FlagPCD=0;par.FlagSSF=0; FlagMinFunc = 0;
options.nLastSteps=0; 
OptimizeBP;  % Perform Optimization
Report.method=method;   reports(i_run)=Report;
options.precond = 1; 
end

if flagSESOP_1
i_run=i_run+1;method='SESOP-1', if FlagConcaveNorm, method=[method '-concave'],end
options.precond = 0;         % 1 - use user-defined  preconditioner, 0 - don't
FlagSESOP=1; par.FlagPCD=0;par.FlagSSF=0; FlagMinFunc = 0;
options.nLastSteps=1; 
OptimizeBP;  % Perform Optimization
Report.method=method;   reports(i_run)=Report;
options.precond = 1; 
end

if flagCG
i_run=i_run+1;method='CG', if FlagConcaveNorm, method=[method '-concave'],end
options.precond = 0; 
FlagSESOP   = 1; par.FlagPCD=0;par.FlagSSF=0;FlagMinFunc = 0;
options.FlagNonlinCG = 1;    % 1 - Polak-Ribiere Nonlinear Conjugate Gradient method; 0 - any other method
options.nLastSteps=1; 
OptimizeBP; % Perform Optimization
options.FlagNonlinCG = 0;    % 1 - Polak-Ribiere Nonlinear Conjugate Gradient method; 0 - any other method
Report.method=method;  reports(i_run)=Report;
options.precond = 1; 
end

if flagSesopLBFGS
i_run=i_run+1;method='SesopLBFGS', if FlagConcaveNorm, method=[method '-concave'],end
options.precond = 0; 
FlagSESOP   = 1; par.FlagPCD=0;par.FlagSSF=0;FlagMinFunc = 0;
options.FlagLBFGS = 1;    % 1 - LBFGS iside sesop code; 0 - any other method
options.nLastSteps=1; 
OptimizeBP; % Perform Optimization
options.FlagLBFGS = 0;    % 1 - Polak-Ribiere Nonlinear Conjugate Gradient method; 0 - any other method
Report.method=method;  reports(i_run)=Report;
options.precond = 1; 
end

if flagSSF_LBFGS
i_run=i_run+1;method='SSF-LBFGS', if FlagConcaveNorm, method=[method '-concave'],end
options.precond = 1; 
FlagSESOP   = 1; par.FlagPCD=0;par.FlagSSF=1;FlagMinFunc = 0;
options.FlagLBFGS = 1;    % 1 - LBFGS iside sesop code; 0 - any other method
options.nLastSteps=1; 
OptimizeBP; % Perform Optimization
options.FlagLBFGS = 0;    % 1 - Polak-Ribiere Nonlinear Conjugate Gradient method; 0 - any other method
Report.method=method;  reports(i_run)=Report;
end

if flagPCD_LBFGS
i_run=i_run+1;method='PCD-LBFGS', if FlagConcaveNorm, method=[method '-concave'],end
options.precond = 1; 
FlagSESOP   = 1; par.FlagPCD=1;par.FlagSSF=0;FlagMinFunc = 0;
options.FlagLBFGS = 1;    % 1 - LBFGS iside sesop code; 0 - any other method
options.nLastSteps=1; 
OptimizeBP; % Perform Optimization
options.FlagLBFGS = 0;    % 1 - Polak-Ribiere Nonlinear Conjugate Gradient method; 0 - any other method
Report.method=method;  reports(i_run)=Report;
end

   
   
if flagPCD
i_run=i_run+1;method='PCD', if FlagConcaveNorm, method=[method '-concave'],end
FlagSESOP   = 1; par.FlagPCD=1;FlagSSF = 0;%options.max_newton_iter = 1;
options.nLastSteps=0; 
OptimizeBP; % Perform Optimization
Report.method=method;  reports(i_run)=Report;
end

if flagPCD_CG
i_run=i_run+1;method='PCD-CG',if FlagConcaveNorm, method=[method '-concave'],end
FlagSESOP   = 1; par.FlagPCD=1;par.FlagSSF = 0;
options.FlagNonlinCG = 1;    % 1 - Polak-Ribiere Nonlinear Conjugate Gradient method; 0 - any other method
options.nLastSteps=1; 
OptimizeBP; % Perform Optimization
options.FlagNonlinCG = 0;    % 1 - Polak-Ribiere Nonlinear Conjugate Gradient method; 0 - any other method
Report.method=method;  reports(i_run)=Report;
end

if flagPCD_SESOP_7
i_run=i_run+1;method='PCD-SESOP-7', if FlagConcaveNorm, method=[method '-concave'],end
FlagSESOP   = 1; par.FlagPCD=1;FlagSSF = 0;
options.nLastSteps=7; 
OptimizeBP; % Perform Optimization
Report.method=method;  reports(i_run)=Report;
end


if flagSSF
i_run=i_run+1;method='SSF', if FlagConcaveNorm, method=[method '-concave'],end
FlagSESOP   = 1; par.FlagPCD=0;par.FlagSSF = 1;
max_newton_iter=options.max_newton_iter;options.max_newton_iter = 0;
options.nLastSteps=0; 
OptimizeBP; % Perform Optimization
Report.method=method;  reports(i_run)=Report;
options.max_newton_iter = max_newton_iter;
end


if flagSSF_lsrch
i_run=i_run+1;method='SSF-lsrch', if FlagConcaveNorm, method=[method '-concave'],end
FlagSESOP   = 1; par.FlagPCD=0;par.FlagSSF = 1;
options.nLastSteps=0; 
OptimizeBP; % Perform Optimization
Report.method=method;  reports(i_run)=Report;
end

if flagSSF_CG
i_run=i_run+1;method='SSF-CG', if FlagConcaveNorm, method=[method '-concave'],end
options.precond = 1; FlagSESOP   = 1; par.FlagPCD=0;par.FlagSSF = 1;
options.FlagNonlinCG = 1;
options.nLastSteps=1; 
OptimizeBP; % Perform Optimization
Report.method=method;  reports(i_run)=Report;
options.FlagNonlinCG = 0;
end

if flagSSF_SESOP_7
i_run=i_run+1;method='SSF-SESOP-7', if FlagConcaveNorm, method=[method '-concave'],end
FlagSESOP   = 1; par.FlagPCD=0;par.FlagSSF = 1;
options.nLastSteps=7; 
OptimizeBP; % Perform Optimization
Report.method=method;  reports(i_run)=Report;
end

if flagLBFGS
i_run=i_run+1;method='L-BFGS', if FlagConcaveNorm, method=[method '-concave'],end
FlagSESOP   = 0; FlagMinFunc = 1;
OptimizeBP; % Perform Optimization
Report.method=method;  reports(i_run)=Report;
FlagMinFunc = 0;
end


if flagL1_LS_IntPoint && ~FlagConcaveNorm
   i_run=i_run+1;method='L1-LS-IntPoint',
   tmp=par.period_show_progress;
   par.period_show_progress=1;
   FlagBoyd=1;FlagSESOP=0;
   OptimizeBP; % Perform Optimization
   Report.method=method;  reports(i_run)=Report;
   FlagBoyd=0;FlagSESOP=1;
   par.period_show_progress=tmp;
end

if flagNesMz_SSF
i_run=i_run+1;method='NesMz-SSF', if FlagConcaveNorm, method=[method '-concave'],end
options.FlagNesMz=1; 
FlagSESOP   = 1; par.FlagPCD=0;par.FlagSSF = 1;
max_newton_iter=options.max_newton_iter;options.max_newton_iter = 0;
options.nLastSteps=0; 
OptimizeBP; % Perform Optimization
Report.method=method;  reports(i_run)=Report;
options.max_newton_iter = max_newton_iter;
options.FlagNesMz=0; 
end

if flagFISTA_SSF
i_run=i_run+1;method='FISTA', if FlagConcaveNorm, method=[method '-concave'],end
options.FlagFISTA=1; 
FlagSESOP   = 1; par.FlagPCD=0;par.FlagSSF = 1;
max_newton_iter=options.max_newton_iter;options.max_newton_iter = 0;
options.nLastSteps=0; 
OptimizeBP; % Perform Optimization
Report.method=method;  reports(i_run)=Report;
options.max_newton_iter = max_newton_iter;
options.FlagFISTA=0; 
end


if flagPRE_CG
i_run=i_run+1;method='PRE-CG',if FlagConcaveNorm, method=[method '-concave'],end
FlagSESOP   = 1; par.FlagPCD=0;FlagSSF = 0;FlagMinFunc = 0;
options.FlagNonlinCG = 1;    % 1 - Polak-Ribiere Nonlinear Conjugate Gradient method; 0 - any other method
options.nLastSteps=1; 
OptimizeBP; % Perform Optimization
options.FlagNonlinCG = 0;    % 1 - Polak-Ribiere Nonlinear Conjugate Gradient method; 0 - any other method
Report.method=method;  reports(i_run)=Report;
end

if flagPRE_SESOP
i_run=i_run+1;method='PRE-SESOP', if FlagConcaveNorm, method=[method '-concave'],end
options.precond = 1;         % 1 - use user-defined  preconditioner, 0 - don't
FlagSESOP=1; par.FlagPCD=0;FlagSSF = 0;FlagMinFunc = 0;
options.nLastSteps=1; 
OptimizeBP;  % Perform Optimization
Report.method=method;   reports(i_run)=Report;
end










%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%              Show results 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ColorPlots=1;

if ColorPlots
   colorvec = [{'k'},{'k--'},{'b'},{'b--'},{'r'},{'r--'},{'g'},{'g--'},{'c'},{'c--'},{'m'},{'m--'},{'k-.'},{'b-.'},{'r-.'},{'g-.'},{'c-.'},{'m-.'}];
else
	colorvec = [{'k'},{'k:'},{'k-.'},{'k--'},{'k.'},{'k+'},{'k^'},{'ko'},{'ks'} {'k*'} {'kv-'}];
end



FlagNewReport=1;

ind=1:i_run;              % which curves to show
report_name=[par.ProblName par.ProblSubName '_' par.ImageName  num2str(n0) '_' sprintf('%.4f',now)];

f_figure_handle=figure('Position',[00 40 400 600],'name',['Objective function' par.ProblName] );

fbest=1e100; for i=1:i_run, fbest=min(fbest,min(reports(i).func_values));end
mylegends=[];
k=0;
for i=ind,
   semilogy(reports(i).nniter_fg, reports(i).func_values-fbest,colorvec{i});
   k=k+1; mylegends{k}=reports(i).method;
    xlabel('Iteration');ylabel('f - fbest');hold on;
end
legend(mylegends);


figure('Position',[400 40 400 600],'name',['Objective with time' report_name] );
for i=ind,
   ttt=reports(i).times;ttt=ttt-ttt(1);
   semilogy(ttt, reports(i).func_values-fbest,colorvec{i});
   xlabel('CPU time, Sec');ylabel('f - fbest'); hold on;
end
legend(mylegends);


figure('Position',[800 40 400 600],'name',['SNR ' report_name] );
for i=ind,
   plot(reports(i).nniter, reports(i).Xsnr,colorvec{i});hold on;
   xlabel('Iteration');ylabel('SNR, Db')
end
legend(mylegends,'Location','SouthEast');grid

figure('Position',[1200 40 400 600],'name',['SNR with time ' report_name] );
for i=1:i_run,
   if isfield(reports(i),'SNRtime') , ttt=reports(i).SNRtime;ttt=ttt-ttt(1);
   else ttt=reports(i).times(reports(i).nniter+1);ttt=ttt-ttt(1);
   end
   plot(ttt, reports(i).Xsnr,colorvec{i});
   xlabel('CPU time, Sec');ylabel('SNR, Db')
   hold on;
end
legend(mylegends,'Location','SouthEast');grid


Y=reshape(par.y,par.proj_size);
X0=reshape(par.SynthesisFunc(c0,par),m,n);
X=reshape(par.SynthesisFunc(c,par),m,n);


%report_file_name=['reports_' report_name  '.mat'], save(report_file_name, 'reports', 'report_name');
%data_file_name=['data_' report_name  '.mat'], save(data_file_name, 'X00', 'Y', 'X0',X);


%beep;pause(0.3);
%tts('Basis Pursuit test has been finished');
%tts('Pozoveetia Meashu beestroo beestroo !');
%for i=1:8, beep;pause(0.1);end

%subplot(224); plot([c00+1,  c-1]); % add +-1 to the data in order to separate the plots
%title('Basis Pursuit results:\newline Original (blue) and recovered (green) coefs');

%profile viewer
