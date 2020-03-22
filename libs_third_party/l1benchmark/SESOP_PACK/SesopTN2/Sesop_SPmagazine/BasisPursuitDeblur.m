% Basis Pursuit Deblur using unconstrained optimization tool SesopTn
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
%
%  Michael  Zibulevsky 07.07.2008; 17.07.2008; 03.08.2008;
%  10.06.2009; 04.08.2009; 23.08.2009; 05.01.2010
%
% Copyright (c) 2008-2010 All rights reserved. No warranty. Free for academic use

init_BP;  % Init parameters

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Choose the methods you want to run (change 0 to 1)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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
flagFISTA_SSF=0;
flagPRE_CG=0;
flagPRE_SESOP=0;


flagPCD=1;
flagSSF=1;
flagPCD_SESOP_7=1;
flagSSF_SESOP_7=1;
flagFISTA_SSF=1;
flagLBFGS=1;

options.PeriodRestoreAx=8;    % period of SESOP restorations to avoid error accumulation
options.PeriodRestoreAp=4;

options.max_newton_iter = 1;   % Max Newton iterations in subspace optimization



i_run=0;clear reports;  % clear counter and reports structure

INIT_DATA=1;              % 1 - init random data;  0 - use old data from the previous run
options.ContinueOldRun=0; % 1- We can run again starting from  achieved solution with changed smoothing/concavity parameters


n0=64;  m0=n0              % Image size: m0-by-n0 


FlagConcaveNorm=0;         % 1 - use concave norn; 0 - convex one


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Choose test image (uncomment one of four lines below)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


par.ImageName='Phantom';
%par.ImageName='Sparse'; FlagSparseImage=1; image_sparsity=0.05; 
%par.ImageName='Lena'; 
%par.ImageName='Peppers'; 

par.ProblName='Deblur'; FlagDenoise=0;   % when FlagDenoise=1 - use delta-function blur (just denoising problem)

par.weight_abs_penalty = 1e-3;  % Weight of smooth abs penalty (mu)
par.eps_smooth_abs     = 1e-3;  % Smoothing parameter of asbsolute value approximation


par.sigma_noise = 1e-2;         % STD of white Gaussian noise in simulated measurements

options.max_sesop_iter  = 100;  % Max  SESOP iterations

options.period_show_progress=10; % Periodicity of showing sesop and user plots





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%        Define Analysis/Synthesis Operator
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if FlagSparseImage,  % For sparse image use identity operator
   par.SynthesisFunc = @(c,par) c;
   par.AnalysisFunc  = @eye_analysis;
else
   par.wavelet_cqf = daubcqf(4,'min');  % Wavelet definitions
   par.wavelet_levels = 4;              %
   par.SynthesisFunc = @TI2Dsynthesismz;% 
   par.AnalysisFunc  = @TI2Danalysismz; % Translation-invariant wavelet transform
end


%options.TESTGRAD=1;  % 1 - Test gradient numerically, plot testing results (without run optimization); 0 - don't test


if INIT_DATA 
   
   create_test_image; % create/read image X00(m0,n0)
   [m,n]=size(X00);
   par.imagesize=[m,n];

   
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %
   %              Prepare  blur kernel
   %
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   if strcmp(par.ProblName,'Deblur')
      %h=gausswin(9); h2=conv2(h,h');     
      [ii,jj]=meshgrid(-7:1:7); h2=1./(ii.^2+jj.^2+1);
      if FlagDenoise, h2=1; end % Delta-function blur (pure denoising problem)
      h2=h2/sum(h2(:));
      par.blur_kernel=h2; %figure;mesh(h2);title('Blur kernel');
      par.proj_size=par.imagesize;
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%       Define Projector Operator and its adjoint
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(par.ProblName,'Deblur')
   par.projector    =@blurring;
   par.projector_adj=@blurring_adj;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Overall Projector-Synthesis Operator and its adjoint
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

par.multA=@projector_synthesis;
par.multAt=@analysis_projector_adj;



if INIT_DATA    &&  ~options.ContinueOldRun ,  
   
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %
   %          Simulate noisy measurements
   %
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   
   [coefs,coefs_size]=feval(par.AnalysisFunc, X00(:),par);
   %[coefs,coefs_size]=par.AnalysisFunc(X00(:),par);      % this line does not work in matlab2004
   par.coefs_size=coefs_size;
   
	par.x00=X00(:);
   par.y00=par.projector(X00(:),par);
   par.y=par.y00+par.sigma_noise*randn(size(par.y00));    %  add noise 
   
   
   
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %
   %   Compute diag_AtA and   max singular value of A
   %
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
   
	fprintf('\n Compute diag_AtA  to be  used for preconditioning...')
	par.diag_AtA=StochasticCalcDiagAtA(par.multAt,size(par.y),100,par); % For large matrices
	%par.diag_AtA=diag(A'*A);                                       % For small matrices
   fprintf('Done\n');
   figure;plot(par.diag_AtA);title('Estimated diagonal of AtA')
   
   c0=par.multAt(par.y,par);
   
   fprintf('\n Compute max singular value of A used in SSF...');
   max_sing_value_A = max_singular_value(par.multA,par.multAt,randn(size(c0)),100,par)
   par.max_sing_value_A=max_sing_value_A;
   fprintf('Done\n');

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Prepare starting point for optimization
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 

c0=par.multAt(par.y,par);
c0=(1/max_sing_value_A^2)*c0;    % Starting point for optimization



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Setup func_u and func_x used in SESOP model
%
%        f(x) = f_u (Ax) + f_x (x)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

par.func_u =@ls_fgh;          % user function  f(u) = 0.5*||u - y||^2
par.func_x=@sum_abs_smooth;   % user function  f(x) = mu*sum(abs_smoothed_eps(x))
%par.func_x=@sum_fast_abs_smooth_new;   % user function  f(x) = mu*sum(abs_smoothed_eps(x))





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%       PERFORM OPTIMIZATION and plot results
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

optimize_all;




