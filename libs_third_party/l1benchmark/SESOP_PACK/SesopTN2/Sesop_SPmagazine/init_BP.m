%init_BP

%addpath('..');
addpath('../');addpath('../util_mz'); addpath('../tomogr_regul');
addpath('../../minFunc')
addpath('../../Rice_Wavelab')
addpath('../../l1_ls_matlab')

global GlobalNiter GlobalGradNorms GlobalFuncValues GlobalTimes GlobalNNiter GlobalNNiterFG GlobalSNRtime GlobalXsnr GlobalY_RNR; 
global Global_conc_lognorm; 
Global_conc_lognorm.simple=1;
Global_conc_lognorm.on = 0;      % 1 - use Lp_norm smooth approximation, p<1, using two log's
%                                  0 - use convex smooth |c| approximation


%maxNumCompThreads(1); % To avoid multiple CPU use (for correct CPU time measuring)
FlagDenoise=0;FlagSparseImage=0; image_sparsity=1;resize_factor=1;

options=sesoptn_optionset;  % Get default options structure (see comments in optionset_sesoptn.m) 

options.PeriodRestoreAx=8;    % period of SESOP restorations to avoid error accumulation
options.PeriodRestoreAp=1;

options.ShowSesopPlots=1;           % 1 - plot iteration progress, 0 - don't

options.dxTol = 1e-16;    % norm of change_in_x
options.dfTol = 1e-16;   % abs change_in_objective_func
options.gTol  = 1e-16;    % norm of the gradient


par.ProblSubName=''; 

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
flagSesopLBFGS_1Nwt=0;
flagSesopLBFGS_Wolfe=0;
flagSSF_LBFGS=0;
flagPCD_LBFGS=0;
flagPCD_LBFGS_1Nwt=0;
flagSesopLBFGS_Wolfe=0;
flagPCD_LBFGS_Wolfe=0;
flagSSF_LBFGS_Wolfe=0;
flagSSF_LBFGS_1Nwt=0;

options.report_func=@deblur_report_func;   % User function reference to display/collect iteration progress;
options.sesop_figure_name=sprintf('SESOPtn'); %d CGstepsPerTN; precond=%d PCD=%d SSF=%d', options.max_iter_CGinTN,options.precond, par.FlagPCD,par.FlagSSF);

par.report_func=options.report_func;