% L1 - L2 Linear Support Vector Machine using SesopTN optimization tool
%
%
%  Let matrix [X,y] contain a data points x_k and +-1 labels y_k at rows k=1,..,m
%  We create matrix  A = [-diag(y)*X, -y];
%
% SVM training via smooth unconstrained optimization:
%
%            min_wb  mu*sum(abs_smooth(w)) + c*svm_penalty(A*[w;b])  
%
% We have a continum between L1 and L2 penalties on separation vector and
% on constraint violation. 
% When smoothing parameters >> 0, the penalties become close to quadratic
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


% Michael Zibulevsky, 05.08.2008        
%
% Copyright (c) 2008. All rights reserved. Free for academic use. No warranty 


%profile on;
addpath('..'); addpath('../util_mz');


global Global_sesop_iter Global_nn_errors; 
Global_sesop_iter=[];Global_nn_errors = [];


INIT_DATA=1;         % 1 - init random data;  0 - use old data from the previous run
CONTINUE_OLD_RUN=0;  % 1- We can run again starting from  achieved solution with changed smoothing/concavity parameters

n=2000;  % data dimension 
m=300; % number of training examples

par.eps_smooth_abs=0.1;          % Smoothing parameter of absolute value of w 
par.eps_smooth_const_linear=0.1; % Smoothing parameter of const-linear penalty on separating band violation
c0=10;                    % Weight of penalty for violation of the separating band
 

options=sesoptn_optionset;  % Get default options structure (see comments in optionset_sesoptn.m) 

options.sesop_figure_name=sprintf('SESOPtn progress,   %d CG steps per TN iteration',options.max_iter_CGinTN);
%options.report_func=@svm_show_progress;
%par.report_figure_handle= figure('Position',[10 40 400 600], 'Name',sprintf('SVM validation error rate'));






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                 Generate data for classification:  
%
%        Two Gaussian clouds shifted by delta in each coodinate
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if INIT_DATA,
	delta=4/sqrt(n);
	
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
	
	
	wb0=randn(n+1,1);  % Initial value of the separation vector and bias
	par.A=DataMatrix;
	par.c=c0*n/m;

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
%                   Perform SESOP-TN optimization
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if CONTINUE_OLD_RUN, wb_init=wb;
else                 wb_init=wb0;
end

par.func_u = @svm_penalty;       % penalty for violation of separating strip
%par.func_u = @svm_penalty_huber; % penalty for violation of separating strip
par.func_x = @svm_diag_penalty; % sum abs_smooth(w)
par.multA  = @(x,par) multMatr1(DataMatrix,x,0);  % user function  y=Ax
par.multAt = @(x,par) multMatr1(DataMatrix,x,1);  % user function  y=A'*x

% 
% par.func_u = @svm_penalty;      % penalty for violation of separating strip
% par.func_x = @svm_diag_penalty; % sum abs_smooth(w)
% par.multA  = @(x,par) multMatr(DataMatrix,x);     % user function  y=Ax
% par.multAt = @(x,par) multMatrAdj(DataMatrix,x);  % user function  y=A'*x


tic
[wb,report]=sesoptn(wb_init(:),par.func_u, par.func_x, par.multA, par.multAt,options,par);
toc

b=wb(end);
w=wb(1:end-1);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% %
% %  Optimization using fminunc.m from Matlab Optimization toolbox 
% %
% %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% disp('Optimization using fminunc.m from Matlab Optimization toolbox')
% options1 =optimset('LargeScale','off','GradObj','on','Display','iter','MaxIter',200);
% c1 =fminunc( @fg, wb_init(:), options1, par);




% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% %
% %  Optimization using minFunc of Mark Schmidt 
% %                              http://www.cs.ubc.ca/~schmidtm
% %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% disp('Optimization using minFunc of Mark Schmidt')
% addpath('../../minFunc')
% options1.Display='final';
% %options1.Display='iter';
% options1.Method = 'lbfgs';
% options1.MaxIter=2000;
% options1.MaxFunEvals=8000;
% options1.TolFun=1e-12;
% options1.TolX=1e-10;
% options1.Corr=8; % Number of previous steps used in L-BFGS update
% 
% tic;
% par.flagXnew=1;
% [wb,f] =minFunc( @fg, wb_init(:), options1,par);
% toc
% fprintf('f_opt=%g ',f)


% %%%%%%% Some testing
% u=DataMatrix*wb;
% %u=rand(size(u));
% z=rand(size(u));
% par.flagXnew=1;
% [f,g,Hz]=svm_penalty(u,z,par)
% [f1,g1,H1z]=svm_penalty_old(u,z,par)
% f-f1,g-g1,Hz-H1z

%profile viewer
