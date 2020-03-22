% OptimizeBP


if  ~options.ContinueOldRun 
   GlobalNNiter=[];GlobalNiter=0;GlobalGradNorms=[]; GlobalFuncValues=[]; GlobalTimes=[];GlobalSNRtime=[];GlobalXsnr=[];GlobalY_RNR=[];
   GlobalNNiterFG=[];
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   User preconditioning function: d_new=mult_precond (-g, x, par)
%   (active, when options.precond=1)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if par.FlagPCD || par.FlagSSF   % PCD or SSF direction 
   options.precond = 1;
	options.mult_precond =@mult_precond_pcd_ssf;

else         % Diagonal preconditioning
	
	options.mult_precond = @mult_diag_precond;
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%      Perform SESOP optimization
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




if options.ContinueOldRun,   c_init=c;
else                   c_init=c0;
end

if FlagSESOP,
   if ~options.ContinueOldRun
      par.BPblur_report_figure_handle=figure('Position',[10 40 600 700],'name',['SESOP ' par.ProblName]);
   end
   tic
   [c,Report]=sesoptn(c_init,par.func_u, par.func_x, par.multA, par.multAt,options,par);
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
options1.Display='final';
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
   if ~options.ContinueOldRun
      par.BPblur_report_figure_handle=figure('Position',[10 40 600 900],'name',['L-BFGS ' par.ProblName] );
   end
   tic;
   [c,f] =minFunc_with_report( @fg, c_init, options1,par);
   toc
   fprintf('f_opt=%g ',f)

   Report.gradnorms   = GlobalGradNorms;
   Report.func_values = GlobalFuncValues;
   Report.times       = GlobalTimes;
   %Report.nniter_fg   = GlobalNNiter;
   

end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%     Interior Point Solver (Boyd)
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if FlagBoyd
   if ~options.ContinueOldRun
      par.BPblur_report_figure_handle=figure('Position',[10 40 600 900],'name',['Boyd_IntPoint ' par.ProblName] );
   end
  
   par.eval_fg=@fg_with_report;
   A = operMultA(par);
   At=A';
   rel_tol = 0.01; % relative target duality gap
   eta   = 1e-3; %scalar; parameter for PCG termination (default: 1e-3)
   quiet=false;   %suppress printing message when true (default: false)

   pcgmaxi =200; % scalar; number of maximum PCG iterations (default: 5000)
   totalpcgmaxi = options.max_sesop_iter;% MZIB: total number of maximum PCG iterations 


   %run the l1-regularized least squares solver
   %profile on
   profile off
   [c,status]=l1_ls_mz(A,At,numel(par.y),numel(c_init),par.y,par.weight_abs_penalty,...
      rel_tol,quiet,eta,pcgmaxi,totalpcgmaxi);
   %profile viewer
   
   Report.gradnorms   = GlobalGradNorms;
   Report.func_values = GlobalFuncValues;
   Report.times       = GlobalTimes;
   %Report.nniter_fg   = GlobalNNiter;
   

end

if isempty(GlobalNNiterFG), Report.nniter_fg=[1:length(Report.func_values)];
else                        Report.nniter_fg=GlobalNNiterFG;
end

Report.nniter      = GlobalNNiter;
Report.Xsnr        = GlobalXsnr;
Report.Yrnr        = GlobalY_RNR;
Report.SNRtime     =GlobalSNRtime;

Report.ProblName=par.ProblName;
Report.ImageName=par.ImageName;
Report.sigma_noise =par.sigma_noise ;
Report.weight_abs_penalty=par.weight_abs_penalty;
Report.eps_smooth_abs=par.eps_smooth_abs;

Report.FlagDenoise=FlagDenoise;
Report.image_sparsity=image_sparsity;
Report.image_size=size(X00);
Report.resize_factor=resize_factor;

Report.options.PeriodRestoreAx=options.PeriodRestoreAx;
Report.options.PeriodRestoreAp=options.PeriodRestoreAp;

Report.FlagNewReport=1;




%GlobalNiter GlobalGradNorms GlobalFuncValues GlobalTimes GlobalNNiter GlobalXsnr GlobalY_RNR
