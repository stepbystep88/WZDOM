
% Michael Zibulevsky, 05.08.2008; 04.02.2009        

addpath('../');addpath('../util_mz');

global Global_sesop_iter; Global_sesop_iter=[];


%profile on;
%profile off;
% 
% Test_problem_number=4;n=2; m=n;  % Rosenbrock
% Test_problem_number=5;n=3; m=n;  % Helical Valley
% Test_problem_number=6;n=4; m=n;  % Powell singular func
% Test_problem_number=7;n=2; m=n;  % Freudenstein & Roth
% 
% 
%Test_problem_number=8;n=3; m=15;   % Bard
% Test_problem_number=15;n=100; m=n+5; % Chebyquad (any n, m>=n)
%Test_problem_number=22;n=50; m=1; % Exp and squares (any n; m=1) % SESOP-TN is better here!!!
Test_problem_number=22;n=200; m=1; % Exp and squares (any n; m=1) % SESOP-TN is better here!!!




 



FlagMinFunc=0; % 1 - use munFunc L-BFGS optimization; 0 - don't



CONTINUE_OLD_RUN=0;  % 1- We can run again starting from  achieved solution with changed smoothing/concavity parameters
%maxCG=[0 10  100];
%PureTN=[0,1];
%maxCG=[10 100];PureTN=[0 1]; max_iter_CG_total=1000;
%maxCG=[10];PureTN=[0]; max_iter_CG_total=200;
%maxCG=[5 10 20 50 100];PureTN=[0 1]; max_iter_CG_total=1500;  
%maxCG=[5 10 20];PureTN=[0 1]; max_iter_CG_total=1500;  
%maxCG=[10 14 16 18 20 100];PureTN=[0 ]; max_iter_CG_total=1000;  
%maxCG=[10 20 50 70 90 100];PureTN=[0 ]; max_iter_CG_total=1000;  
%maxCG=[10 40 100];PureTN=[0 1]; max_iter_CG_total=3000;  
maxCG=[1 10 40];PureTN=[0 1]; max_iter_CG_total=3000;  


 

options=sesoptn_optionset;  % Get default options structure (see comments in optionset_sesoptn.m) 

options.AnalyticHessXmult=1;   % 1 Analytic multiplication by Hess in func_x; 0 - by finite differences

options.PureTN=0;             % 1 - Pure Truncated Newton; 0 - Truncated Newton with  SESOP (recommended); 

options.FlagNonlinCG    = 0;  % 1 - Polak-Ribiere Nonlinear Conjugate Gradient method; 0 - any other method


options.max_newton_iter=5;     % Max Newton iterations in subspace optimization
options.max_iter_CGinTN = 10;  % Conj.Grad steps in Truncated Newton (if  set to 0, TN  not used)
options.max_sesop_iter=max_iter_CG_total/(options.max_iter_CGinTN+1);    % Max  SESOP iterations


options.nLastSteps=1;         % SESOP subspace will include: nLastSteps + [precond] gradient + [3 TN directions] + [2 Nemirovski directions])
options.FlagGrad_dir=1;       % 1 include [precond] gradient into subspace; 0 - don't
options.FlagNemirovski_dir=0; % 1 - Add two directions of Nemirovski: x-x0 and  sum w_k*g_k

options.PeriodRestoreAx=1;
options.PeriodRestoreAp=1;

options.dxTol=-1e-60;         % stopping criteria norm(change_in_x)
options.dfTol=-1e-60;         % abs(change_in_objective_func)
options.gTol=1e-6;            % norm(gradient)

																
options.ShowSesopPlots=1;             % 1 - plot iteration progress, 0 - don't
options.period_show_progress=200;     % Periodicity of showing sesop and user plots
options.sesop_figure_name=sprintf('SESOPtn progress,   %d CG steps per TN iteration',options.max_iter_CGinTN);
%options.report_func=@svm_show_progress;
%par.report_figure_handle= figure('Position',[10 40 400 600], 'Name',sprintf('SVM validation error rate'));


par.flagXnew=1; % Needed for fg.m used by minFunc

ColorPlots=0;

if ColorPlots
	%colorvec = [{'k'},{'b'},{'r'},{'g'},{'y'},{'c'},{'m'},{'.'},{'--'},{'r.'}];
	colorvec = [{'k'},{'k:'},{'b'},{'b:'},{'r'},{'r:'},{'g'},{'g:'},{'y'},{'c'},{'m'},{'.'},{'--'},{'r.'}];
	
else
	%colorvec = [{'k:'},{'k-.'},{'k--'},{'k.'},{'k'},{'k+'},{'k^'},{'ks'},{'k'}];
	%colorvec = [{'k'},{'k:'},{'k-.'},{'k--'},{'k.'},{'k+'},{'k^'},{'ko'},{'ks'}];
	%colorvec = [{'k'},{'k:'},{'k-.'},{'ko'},{'k+'},{'k^'},{'ko'},{'ks'},{'k--'}];
	colorvec = [{'k'},{'k:'},{'k.'},{'k+'},{'k-.'},{'k^'},{'ko'},{'ks'},{'k--'}];
	%colorvec = [{'k'},{'k.'},{'ko'},{'k*'},{'k^'},{'ks'},{'kd'},{'k>'},{'kv'}];
end







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%                   Perform SESOP-TN optimization
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[par, x_init] = uctpget(Test_problem_number,m,n);
par.pt=1;


% tic
% [x,report]=sesoptn_svm(x_init(:),[], @uctpval_sesop, [],[],options,par);
% toc

%profile viewer

 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%  Optimization using minFunc of Mark Schmidt 
%                              http://www.cs.ubc.ca/~schmidtm
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if FlagMinFunc
	%profile on;

	disp('Optimization using minFunc of Mark Schmidt')
	addpath('../../../minFunc')
	options1 =optimset('GradObj','on','Display','iter','MaxIter',800);
	options1.Display='final';
	options1.Method = 'lbfgs';
	options1.TolFun=1e-10;
	options1.TolX=1e-10;
	options1.Corr=8; % Number of previous steps used in L-BFGS update

	tic;
	[x,f] =minFunc( @uctpval, x_init(:), options1,par);
	toc
	fprintf('f_opt=%g ',f)
	%profile viewer
end

 
  
 
 

kkk=length(maxCG);


j=0;
for i=1:kkk,
	for k=PureTN
		options.PureTN=k;
		j=j+1;
		fprintf('\n max_iter_CGinTN=%d',maxCG(i));
		options.max_iter_CGinTN=maxCG(i)+options.PureTN;
		options.max_sesop_iter=max_iter_CG_total/(options.max_iter_CGinTN+1);
		options.sesop_figure_name=sprintf('SESOPtn progress,   %d CG steps per TN iteration',options.max_iter_CGinTN);
		%options.period_show_progress=max(1,ceil(1/(maxCG(i)+1)));     % Periodicity of showing sesop and user plots
		[x,report(j)]=sesoptn(x_init(:),[], @uctpval_sesop, [],[],options,par);
	end
end


  

fbest=1e100;
for i=1:kkk,fbest=min(fbest,min(report(i).func_values));end

fighandle1=figure('name', sprintf('Func values with iterations; Fbest=%.8g',fbest)); 
fighandle2=figure('name', sprintf('Function values with time, PureTN=%d',options.PureTN)); 
legendstrings=[];

j=0;
for i=1:kkk,
	for k=PureTN
		options.PureTN=k;
		j=j+1;
		ff=report(j).func_values;
		nncg=[0:length(ff)-1]*(maxCG(i) +1);
		figure(fighandle1);
		%semilogy(nncg,ff,colorvec{i}); hold on
		%subplot(121);
		semilogy(nncg,ff-fbest,colorvec{j}); hold on
		%subplot(122);semilogy(nncg,report(j).gradnorms,colorvec{j}); hold on
		%legendstrings{j}=sprintf('%d CG steps in TN',maxCG(i));

		figure(fighandle2);
		subplot(121);semilogy(report(j).times,ff-fbest,colorvec{j}); hold on
		subplot(122);semilogy(report(j).times,report(j).gradnorms,colorvec{j}); hold on
	end
 	legendstrings{j-1}=sprintf('SESOP-TN:   %d CG steps',maxCG(i));
 	legendstrings{j}  =sprintf('Classic TN: %d CG steps',maxCG(i));

end
%figure(fighandle1);subplot(121);legend(legendstrings);%subplot(122);legend(legendstrings);
%figure(fighandle2);subplot(121);legend(legendstrings);%subplot(122);legend(legendstrings);
figure(fighandle1);legend(legendstrings,'Location','SouthWest');
xlabel('Global CG iteration count'); ylabel('Objective function residual')

figure(fighandle2);legend(legendstrings);






























% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% %
% %  Optimization using fminunc.m from Matlab Optimization toolbox 
% %
% %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% disp('Optimization using fminunc.m from Matlab Optimization toolbox')
%options1 =optimset('LargeScale','off','GradObj','on','Display','iter','MaxIter',200);
% c1 =fminunc( @uctpval, x_init(:), options1, par);






%for i=1:5, beep;pause(0.01*i);end