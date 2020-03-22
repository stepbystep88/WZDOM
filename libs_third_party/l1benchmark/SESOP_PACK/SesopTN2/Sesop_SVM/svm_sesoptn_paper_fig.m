% L1 - L2 Linear Support Vector Machine using SesopTN optimization tool
%
%
%  Let matrix [X,y] contain data points x_k and +-1 labels y_k at rows k=1,..,m
%  We create matrix  A = [-diag(y)*X, -y];
%
% SVM training via smooth unconstrained optimization:
%
%            min_wb  mu*sum(abs_smooth(w)) + c*svm_penalty(A*[w;b])  
%
% We have a continum between L1 and L2 penalties on separation vector and
% on constraint violation. 
% When smoothing parameters >> 1, the penalties become close to quadratic
%
%par.eps_smooth_abs          - Smoothing parameter of absolute value of w 
%par.eps_smooth_const_linear - Smoothing parameter of const-linear penalty on separating band violation
%par.c                       - Weight of penalty for violation of the separating band
%
%
% See also:
% Narkiss, G. and Zibulevsky, M. (2005). "Support Vector Machine via Sequential 
% Subspace Optimization", Tech. Report CCIT No 557, EE Dept., Technion.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% Michael Zibulevsky, 05.08.2008    03.02.2009    
%
% Copyright (c) 2008. All rights reserved. Free for academic use. No warranty 


%profile on;
%profile off;


addpath('../');addpath('../util_mz');

%addpath(genpath('../../light_svm_mex601'));
%addpath('../../mex_svm_perf-w32')


global Global_sesop_iter Global_nn_errors; 
Global_sesop_iter=[];Global_nn_errors = [];

FlagMinFunc=0; % 1 - use munFunc L-BFGS optimization; 0 - don't


INIT_DATA=1;         % 1 - init data;  0 - use old data from the previous run
READ_DATA =0;        % 1 - read SVM data file 0 - generate Gaussian clouds 
train_part=1/2;      % Part of the data used for training

CONTINUE_OLD_RUN=0;  % 1- We can run again starting from  achieved solution with changed smoothing/concavity parameters
%maxCG=[0 10  100];
%PureTN=[0,1];
maxCG=[1 10 40];PureTN=[0 1]; max_iter_CG_total=130;
%maxCG=[10];PureTN=[0]; max_iter_CG_total=200;
%maxCG=[0];PureTN=[0]; max_iter_CG_total=70;

%data_dir_name='../../../SVM_Demo/svm_data/';
%data_file_name='australian_scale.htm';   FlagSparseFileFormat=1;
%data_file_name='dna.scale';              FlagSparseFileFormat=1;


data_dir_name='../../../../svm_Thorsten_data/';
%data_file_name='CCAT_23149.dat';              FlagSparseFileFormat=1;
data_file_name='astro-ph_29882.dat';          FlagSparseFileFormat=1; train_part=1/20;

%data_dir_name='../../../SVM_Demo/data/';
%data_file_name='phy_train.dat';          FlagSparseFileFormat=0;
%data_file_name='covtype.data';          FlagSparseFileFormat=-1;

n=50;  % data dimension               (if Gaussians to be generated)
m=100;  % number of training examples (if Gaussians to be generated)


par.eps_smooth_abs=5;          % Smoothing parameter of absolute value of w 
par.eps_smooth_const_linear=0.5;  % Smoothing parameter of const-linear penalty on separating band violation
c0=10                           % See below par.c=c0*n/mtrain, Weight of penalty for violation of the separating band
 

options=sesoptn_optionset;  % Get default options structure (see comments in optionset_sesoptn.m) 

options.PureTN=0;             % 1 - Pure Truncated Newton; 0 - Truncated Newton with  SESOP (recommended); 
options.FlagRestrictedMatrix_TN=1; % 1 - restrict TN matrix; 0 - don't

options.FlagNonlinCG    = 0;  % 1 - Polak-Ribiere Nonlinear Conjugate Gradient method; 0 - any other method


options.max_newton_iter=5;     % Max Newton iterations in subspace optimization
options.max_iter_CGinTN = 20;  % Conj.Grad steps in Truncated Newton (if  set to 0, TN  not used)
options.max_sesop_iter=1000/(options.max_iter_CGinTN+1);    % Max  SESOP iterations


options.nLastSteps=2;         % SESOP subspace will include: nLastSteps + [precond] gradient + [3 TN directions] + [2 Nemirovski directions])
options.FlagGrad_dir=1;       % 1 include [precond] gradient into subspace; 0 - don't
options.FlagNemirovski_dir=0; % 1 - Add two directions of Nemirovski: x-x0 and  sum w_k*g_k

options.PeriodRestoreAx=1;
options.PeriodRestoreAp=1;

options.dxTol=-1e-60;         % stopping criteria norm(change_in_x)
options.dfTol=-1e-60;         % abs(change_in_objective_func)
options.gTol=1e-6;            % norm(gradient)

																
options.ShowSesopPlots=1;           % 1 - plot iteration progress, 0 - don't
options.period_show_progress=1;     % Periodicity of showing sesop and user plots
options.sesop_figure_name=sprintf('SESOPtn progress,   %d CG steps per TN iteration',options.max_iter_CGinTN);
options.report_func=@svm_show_progress;
par.report_figure_handle= figure('Position',[10 40 400 600], 'Name',sprintf('SVM validation error rate'));

par.flagXnew=1; % Needed for fg.m used by minFunc

ColorPlots=1;

if ColorPlots
	%colorvec = [{'k'},{'b'},{'r'},{'g'},{'y'},{'c'},{'m'},{'.'},{'--'},{'r.'}];
	colorvec = [{'k'},{'k:'},{'b'},{'b:'},{'r'},{'r:'},{'g'},{'g:'},{'y'},{'c'},{'m'},{'.'},{'--'},{'r.'}];
	
else
	%colorvec = [{'k:'},{'k-.'},{'k--'},{'k.'},{'k'},{'k+'},{'k^'},{'ks'},{'k'}];
	%colorvec = [{'k'},{'k:'},{'k-.'},{'k--'},{'k.'},{'k+'},{'k^'},{'ko'},{'ks'}];
	%colorvec = [{'k'},{'k:'},{'k.'},{'k+'},{'k-.'},{'k^'},{'ko'},{'ks'},{'k--'}];
	colorvec = [{'k'},{'k.'},{'k:'},{'k+'},{'k-.'},{'k^'},{'ko'},{'ks'},{'k--'}];

	%colorvec = [{'k'},{'k:'},{'k-.'},{'ko'},{'k+'},{'k^'},{'ko'},{'ks'},{'k--'}];
	%colorvec = [{'k'},{'k.'},{'ko'},{'k*'},{'k^'},{'ks'},{'kd'},{'k>'},{'kv'}];
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                 Generate data for classification:  
%
%        Two Gaussian clouds shifted by delta in each coodinate
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



if INIT_DATA,
	
	if READ_DATA,
		disp('Reading data...')
		if FlagSparseFileFormat == 1
			[y0,X]=read_sparse([data_dir_name data_file_name]);
		else
			tmp=load('-ascii', [data_dir_name data_file_name]);
			if FlagSparseFileFormat == 0
				y0=tmp(:,2);
				X=tmp(:,3:end);
			else
				y0=tmp(:,end);
				X=tmp(:,2:end-1);
			end
			clear tmp;
		end
		scale_x=max(abs(X))+1e-40;
		X= mulmd(X,1./scale_x); % Scaling to range [-1 1]
		
		disp('Preparing training and validation matrices...')
		y=sign(y0-mean(y0));  % In case of not (+-1) values, force them to +-1
		[M,n]=size(X);
		ind=randperm(M);
		mtrain=ceil(M*train_part); mtest=M-mtrain;
		Xtrain=X(ind(1:mtrain),:);
		ytrain=y(ind(1:mtrain));
		DataMatrix= [ -muldm(ytrain,Xtrain), -ytrain];
		par.X_validation=X(ind(1+mtrain:end),:)';
		par.y_validation=y(ind(1+mtrain:end));
		
		clear X Xtrain;
		
	else  % Generate two Gaussian clouds

		delta=4/sqrt(n);  % Clouds separation shift

		%%%%%%%%%%% Training Set %%%%%%%%%%%%%%%%%%
		m1=floor(m/4);  X1=randn(n,m1); y1=ones(m1,1);
		m2=m-m1; X2=randn(n,m2) + delta*ones(n,m2); y2=-ones(m2,1);
		X=[X1 X2]; y=[y1;y2];
		ind=randperm(m); X=X(:,ind);y=y(ind);
		DataMatrix= [ -mulmd(X,y)', -y];

		%%%%%%%%% Validation Set %%%%%%%%%%%%%%%%%%%
		X1=randn(n,m/2); y1=ones(m/2,1);
		X2=randn(n,m/2) + delta*ones(n,m/2); y2=-ones(m/2,1);
		par.X_validation=[X1 X2]; par.y_validation=[y1;y2];
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	end

	par.c=c0*n/mtrain;
	par.A=DataMatrix;
	wb0=randn(n+1,1);  % Initial value of the separation vector and bias

end



%%%%%%%%%%%%  For debugging
%ind1=[5: m];                                   
%DataMatrix= [mulmd(X,y)', y];
%DataMatrix= [-mulmd(X(:,ind1),y(ind1))', -y(ind1)];
%DataMatrix=DataMatrix(ind1,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%figure;imagesc(x00);colorbar;return



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%                   SVM_light training
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if 0,
	model = svmlearn(X,y,'-t 2 -g 0.3 -c 0.5')
	model=svmperflearn(X,y,'-c 1 -w 3 --t 2 --b 0');
	% with this, you can get information of the model obtained

	struct(model)
	% and if you want to access the support vectors:
	model.supvec
	% or the alphas (multiplied by the class)
	model.alphas
	% if you want to classify more data with the model obtained previously
	predictions=svmperfclassify(more_data,more_labels,model,'-v 1');
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%                   Perform SESOP-TN optimization
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if CONTINUE_OLD_RUN, wb_init=wb;
else                 wb_init=wb0;
end

%par.func_u = @svm_penalty;       % penalty for violation of separating strip
par.func_u = @svm_penalty_huber; % penalty for violation of separating strip
par.func_x = @svm_diag_penalty; % sum abs_smooth(w)
par.multA  = @(x,par) multMatr1(DataMatrix,x,0);  % user function  y=Ax
par.multAt = @(x,par) multMatr1(DataMatrix,x,1);  % user function  y=A'*x

% tic
%[wb,report]=sesoptn(wb_init(:),par.func_u, par.func_x, par.multA, par.multAt,options,par);
% toc

% tic
% [wb,report]=sesoptn_svm(wb_init(:),par.func_u, par.func_x, par.multA, par.multAt,options,par);
% %[wb,report]=sesoptn_svm_0309_2008(wb_init(:),par.func_u, par.func_x, par.multA, par.multAt,options,par);
% toc

% b=wb(end);
% w=wb(1:end-1);



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
		options.period_show_progress=max(1,ceil(1/(maxCG(i)+1)));     % Periodicity of showing sesop and user plots

		[wb,report(j)]=sesoptn(wb_init(:),par.func_u, par.func_x, par.multA, par.multAt,options,par);
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




% j=0;
% for i=1:kkk,
% 	for k=PureTN
% 		options.PureTN=k;
% 		j=j+1;
% 		ff=report(j).func_values;
% 		nncg=[0:length(ff)-1]*(maxCG(i) +1);
% 		figure(fighandle1);
% 		%semilogy(nncg,ff,colorvec{i}); hold on
% 		subplot(121);semilogy(nncg,ff-fbest,colorvec{j}); hold on
% 		subplot(122);semilogy(nncg,report(j).gradnorms,colorvec{j}); hold on
% 		legendstrings{j}=sprintf('%d CG steps in TN',maxCG(i));
% 
% 		figure(fighandle2);
% 		subplot(121);semilogy(report(j).times,ff-fbest,colorvec{j}); hold on
% 		subplot(122);semilogy(report(j).times,report(j).gradnorms,colorvec{j}); hold on
% 	end
% end
% %figure(fighandle1);subplot(121);legend(legendstrings);%subplot(122);legend(legendstrings);
% %figure(fighandle2);subplot(121);legend(legendstrings);%subplot(122);legend(legendstrings);
% figure(fighandle1);legend(legendstrings);
% xlabel('Global CG iteration count'); ylabel('Objective Function')
% 
% figure(fighandle2);legend(legendstrings);
% 





























% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% %
% %  Optimization using fminunc.m from Matlab Optimization toolbox 
% %
% %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% disp('Optimization using fminunc.m from Matlab Optimization toolbox')
%options1 =optimset('LargeScale','off','GradObj','on','Display','iter','MaxIter',200);
% c1 =fminunc( @fg, wb_init(:), options1, par);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%  Optimization using minFunc of Mark Schmidt 
%                              http://www.cs.ubc.ca/~schmidtm
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if FlagMinFunc
	profile on;
	par.multA  = @(x,par) multMatr1(DataMatrix,x,0);  % user function  y=Ax
	par.multAt = @(x,par) multMatr1(DataMatrix,x,1);  % user function  y=A'*x

	disp('Optimization using minFunc of Mark Schmidt')
	addpath('../../minFunc')
	options1 =optimset('GradObj','on','Display','iter','MaxIter',800);
	%options1.Display='final';
	options1.Method = 'lbfgs';
	options1.TolFun=1e-10;
	options1.TolX=1e-10;
	options1.Corr=8; % Number of previous steps used in L-BFGS update

	tic;
	[wb,f] =minFunc( @fg, wb_init(:), options1,par);
	toc
	fprintf('f_opt=%g ',f)
	profile viewer
end

% %%%%%%% Some testing
% u=DataMatrix*wb;
% %u=rand(size(u));
% z=rand(size(u));
% par.flagXnew=1;
% [f,g,Hz]=svm_penalty(u,z,par)
% [f1,g1,H1z]=svm_penalty_old(u,z,par)
% f-f1,g-g1,Hz-H1z

u=par.multA(wb,par); 
v= (u>-1 & u< -1+par.eps_smooth_const_linear); 
sparsity_v=nnz(v)/length(u)
%figure;plot(v)

%profile viewer
for i=1:5, beep;pause(0.01*i);end