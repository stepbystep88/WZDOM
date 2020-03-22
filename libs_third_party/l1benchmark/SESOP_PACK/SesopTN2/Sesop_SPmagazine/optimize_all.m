% optimize_all.m


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%
%       PERFORM Basis Pursuit OPTIMIZATION: calls OptimizeBP script 
%
%        with different methods repeatedly and plots all results
%
%               (called from BasisPursuitDeblur, ...)
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



FlagMinFunc = 0;     % Optimization using L-BFGS, minFunc of Mark Schmidt
FlagSESOP   = 1;     % Optimization using SESOP
FlagBoyd    = 0;     % Optimization using Boyd's l1_ls interior point solver

options.precond = 1;         % 1 - use user-defined  preconditioner, 0 - don't
par.FlagPCD     = 0;
par.FlagSSF     = 0;
par.FlagPCDSSF_inside_TN=0;
options.FlagNonlinCG = 0;    % 1 - Polak-Ribiere Nonlinear Conjugate Gradient method; 0 - any other method
options.max_iter_CGinTN = 0; 

if INIT_DATA, options.ContinueOldRun=0;end
par.period_show_progress=options.period_show_progress;



%i_run=0;clear reports;  % clear counter and reports structure




if FlagConcaveNorm
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %
   %              Case of  concave penalty
   %
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   
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
   c_conc=3;    % argument scaling
   b_conc=100;  %1./par.eps_smooth_abs; % second derivative at the origin
   [p_tmp, h_tmp, b_tmp] = concave_lognorm_preset(b_conc,  c_conc);

   c0=c;  % Use "warm start" given by convex optimization

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
end

if flagSESOP_1
i_run=i_run+1;method='SESOP-1', if FlagConcaveNorm, method=[method '-concave'],end
options.precond = 0;         % 1 - use user-defined  preconditioner, 0 - don't
FlagSESOP=1; par.FlagPCD=0;par.FlagSSF=0; FlagMinFunc = 0;
options.nLastSteps=1; 
OptimizeBP;  % Perform Optimization
Report.method=method;   reports(i_run)=Report;
end

if flagCG
i_run=i_run+1;method='CG', if FlagConcaveNorm, method=[method '-concave'],end
options.precond = 0; 
FlagSESOP   = 1; par.FlagPCD=0;par.FlagSSF=0;FlagMinFunc = 0;
options.FlagNonlinCG = 1;    % 1 - Polak-Ribiere Nonlinear Conjugate Gradient method; 0 - any other method
options.nLastSteps=1; 
OptimizeBP; % Perform Optimization
options.FlagNonlinCG = 0;   
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


if flagFISTA_SSF
i_run=i_run+1;method='FISTA', if FlagConcaveNorm, method=[method '-concave'],end
options.FlagFISTA=1; 
FlagSESOP   = 1; par.FlagPCD=0;par.FlagSSF = 1;
max_newton_iter=options.max_newton_iter; options.max_newton_iter = 0;
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




% if flagNesMz_SSF
% i_run=i_run+1;method='NesMz-SSF', if FlagConcaveNorm, method=[method '-concave'],end
% options.FlagNesMz=1; 
% FlagSESOP   = 1; par.FlagPCD=0;par.FlagSSF = 1;
% max_newton_iter=options.max_newton_iter;options.max_newton_iter = 0;
% options.nLastSteps=0; 
% OptimizeBP; % Perform Optimization
% Report.method=method;  reports(i_run)=Report;
% options.max_newton_iter = max_newton_iter;
% options.FlagNesMz=0; 
% end
%

if flagSesopLBFGS
i_run=i_run+1;method='SesopLBFGS', if FlagConcaveNorm, method=[method '-concave'],end
options.precond = 0; 
FlagSESOP   = 1; par.FlagPCD=0;par.FlagSSF=0;FlagMinFunc = 0;
options.FlagLBFGS = 1;   % 1 - LBFGS iside sesop code; 0 - any other method
options.nLastSteps=0; 
OptimizeBP;              % Perform Optimization
options.FlagLBFGS = 0;   
Report.method=method;  reports(i_run)=Report;
end

if flagSesopLBFGS_1Nwt
i_run=i_run+1;method='SesopLBFGS-1Nwt', if FlagConcaveNorm, method=[method '-concave'],end
options.precond = 0; 
max_newton_iter=options.max_newton_iter;options.max_newton_iter = 1;
FlagSESOP   = 1; par.FlagPCD=0;par.FlagSSF=0;FlagMinFunc = 0;
options.FlagLBFGS = 1;  options.FlagWolfeLineSearch=0; % 1 - LBFGS iside sesop code; 0 - any other method
options.nLastSteps=0; 
OptimizeBP;              % Perform Optimization
options.FlagLBFGS = 0;  options.FlagWolfeLineSearch=0; options.max_newton_iter = max_newton_iter; 
Report.method=method;  reports(i_run)=Report;
end

if flagSesopLBFGS_Wolfe
i_run=i_run+1;method='SesopLBFGS-Wolfe', if FlagConcaveNorm, method=[method '-concave'],end
options.precond = 0; 
max_newton_iter=options.max_newton_iter;options.max_newton_iter = 0;
FlagSESOP   = 1; par.FlagPCD=0;par.FlagSSF=0;FlagMinFunc = 0;
options.FlagLBFGS = 1;  options.FlagWolfeLineSearch=1; % 1 - LBFGS iside sesop code; 0 - any other method
options.nLastSteps=0; 
OptimizeBP;              % Perform Optimization
options.FlagLBFGS = 0;  options.FlagWolfeLineSearch=0; options.max_newton_iter = max_newton_iter; 
Report.method=method;  reports(i_run)=Report;
end

if flagPCD_LBFGS
i_run=i_run+1;method='PCD-LBFGS', if FlagConcaveNorm, method=[method '-concave'],end
options.precond = 1; 
FlagSESOP   = 1; par.FlagPCD=1;par.FlagSSF=0;FlagMinFunc = 0;
options.FlagLBFGS = 1;    % 1 - LBFGS iside sesop code; 0 - any other method
options.nLastSteps=0; 
OptimizeBP; % Perform Optimization
options.FlagLBFGS = 0;    % 1 - Polak-Ribiere Nonlinear Conjugate Gradient method; 0 - any other method
Report.method=method;  reports(i_run)=Report;
end


if flagPCD_LBFGS_1Nwt
i_run=i_run+1;method='PCD-LBFGS-1Nwt', if FlagConcaveNorm, method=[method '-concave'],end
options.precond = 1; 
max_newton_iter=options.max_newton_iter;options.max_newton_iter = 1;
FlagSESOP   = 1; par.FlagPCD=1;par.FlagSSF=0;FlagMinFunc = 0;
options.FlagLBFGS = 1;  options.FlagWolfeLineSearch=0; % 1 - LBFGS iside sesop code; 0 - any other method
options.nLastSteps=0; 
OptimizeBP;              % Perform Optimization
options.FlagLBFGS = 0;  options.FlagWolfeLineSearch=0; options.max_newton_iter = max_newton_iter; 
Report.method=method;  reports(i_run)=Report;
end

if flagPCD_LBFGS_Wolfe
i_run=i_run+1;method='PCD-LBFGS-Wolfe', if FlagConcaveNorm, method=[method '-concave'],end
options.precond = 1; 
max_newton_iter=options.max_newton_iter;options.max_newton_iter = 0;
FlagSESOP   = 1; par.FlagPCD=1;par.FlagSSF=0;FlagMinFunc = 0;
options.FlagLBFGS = 1;  options.FlagWolfeLineSearch=1; % 1 - LBFGS iside sesop code; 0 - any other method
options.nLastSteps=0; 
OptimizeBP;              % Perform Optimization
options.FlagLBFGS = 0;  options.FlagWolfeLineSearch=0; options.max_newton_iter = max_newton_iter; 
Report.method=method;  reports(i_run)=Report;
end

if flagSSF_LBFGS
i_run=i_run+1;method='SSF-LBFGS', if FlagConcaveNorm, method=[method '-concave'],end
options.precond = 1; 
FlagSESOP   = 1; par.FlagPCD=0;par.FlagSSF=1;FlagMinFunc = 0;
options.FlagLBFGS = 1;    % 1 - LBFGS iside sesop code; 0 - any other method
options.nLastSteps=0; 
OptimizeBP; % Perform Optimization
options.FlagLBFGS = 0;    % 1 - Polak-Ribiere Nonlinear Conjugate Gradient method; 0 - any other method
Report.method=method;  reports(i_run)=Report;
end



if flagSSF_LBFGS_1Nwt
i_run=i_run+1;method='SSF-LBFGS-1Nwt', if FlagConcaveNorm, method=[method '-concave'],end
options.precond = 1; 
max_newton_iter=options.max_newton_iter;options.max_newton_iter = 1;
FlagSESOP   = 1; par.FlagPCD=0;par.FlagSSF=1;FlagMinFunc = 0;
options.FlagLBFGS = 1;  options.FlagWolfeLineSearch=0; % 1 - LBFGS iside sesop code; 0 - any other method
options.nLastSteps=0; 
OptimizeBP;              % Perform Optimization
options.FlagLBFGS = 0;  options.FlagWolfeLineSearch=0; options.max_newton_iter = max_newton_iter; 
Report.method=method;  reports(i_run)=Report;
end


if flagSSF_LBFGS_Wolfe
i_run=i_run+1;method='SSF-LBFGS-Wolfe', if FlagConcaveNorm, method=[method '-concave'],end
options.precond = 1; 
max_newton_iter=options.max_newton_iter;options.max_newton_iter = 0;
FlagSESOP   = 1; par.FlagPCD=0;par.FlagSSF=1;FlagMinFunc = 0;
options.FlagLBFGS = 1;  options.FlagWolfeLineSearch=1; % 1 - LBFGS iside sesop code; 0 - any other method
options.nLastSteps=0; 
OptimizeBP;              % Perform Optimization
options.FlagLBFGS = 0;  options.FlagWolfeLineSearch=0; options.max_newton_iter = max_newton_iter; 
Report.method=method;  reports(i_run)=Report;
end



   







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%
%
%             Show results 
%
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






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

figure('Position',[1000 40 400 600],'name',['SNR with time ' report_name] );
for i=1:i_run,
   if strcmp(reports(i).method,'L1-LS-IntPoint') %isfield(reports(i),'SNRtime') 
      ttt=reports(i).SNRtime;ttt=ttt-ttt(1);
   else
      ttt=reports(i).times(reports(i).nniter+1);ttt=ttt-ttt(1);
   end
   plot(ttt, reports(i).Xsnr,colorvec{i});
   xlabel('CPU time, Sec');ylabel('SNR, Db')
   hold on;
end
legend(mylegends,'Location','SouthEast');grid


Y=reshape(par.y,par.proj_size);
X0=reshape(par.SynthesisFunc(c0,par),m,n);
X=reshape(par.SynthesisFunc(c,par),m,n);


%tts('Basis Pursuit test has been finished');
%tts('Pozoveetia Meashu beestroo beestroo !');
%for i=1:8, beep;pause(0.1);end
