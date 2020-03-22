
function vcwb = nntrain_mz(Xtrain,ytrain,Xtest,ytest,par)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%           Random initialization of Neutal Net (NN)  weights  v,W,c,b  
%
%                where     NN(X; vcwb)  =  v'phi(WX+b*1') + c 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


global GlobalErrTrain GlobalErrTest GlobalNNiter
GlobalErrTrain=[]; GlobalErrTest =[]; GlobalNNiter=[];

FlagSesop =   1;   %  1 - sesoptn tool;   
FlagFminunc = 0;   %  1 - fminunc BFGS Quasi-Newton (0ptimization toolbox of MATAB)
FlagMinFunc = 0;   %  1 - minFunc L-BFGS Quasi-Newton


options=sesoptn_optionset;        % Get default options structure (see comments in sesoptn_optionset.m) 
options.max_sesop_iter  = 300;    % Max  SESOP iterations
options.period_show_progress=10;  % Periodicity of showing sesop and user plots
options.PeriodRestoreAx=8;        % To avoid error accumulation
options.PeriodRestoreAp=4;
options.nLastSteps=5;             % SESOP subspace
options.max_newton_iter = 1;      % Max Newton iterations in subspace optimization
options.max_iter_CGinTN=0;


par.eps_sigmoid=0.7;        % Sigmoid parameter for sigmoid_mz function
par.quadrpenpar=1e-1;       % Quadratic penalty for NN weights: used for regularization, to avoid overfitting
par.nneurons=5;            % Number of neurons

par.nnreport_figure_handle= figure('Position',[10 40 400 600], 'Name',sprintf('NN error function: penpar=%g',par.quadrpenpar));



if 1 %INIT_VARIABLES   

	[N,K]=size(Xtrain);
	par.Ktrain_samples=K;
	M=par.nneurons;
	
	v0=(1/sqrt(M))*randn(M,1);
	c0=0;
	W0=(1/sqrt(N))*randn(M,N);
	b0=0.1*randn(M,1);
	vcwb0=[v0;c0;W0(:);b0];
	vcwb=vcwb0;
end

vcwb_init=vcwb;

% if continue_old,	vcwb_init=vcwb;      %  Skip weights initialization, if you want just improve previous optimization results
% else                vcwb_init=vcwb0;
% end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                Learning of Neural Nnet  weights via optimization:   
%
%     min_vcwb  1/2 ||Ytrain - NN(Xtrain;vcwb)||^2 +  1/2 par.quadrpenpar*||vcwb||^2
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




fprintf('\n Perform learning of Neural Net  weights via optimization:\n \n ');

	par.func_x=@diag_quadr_penalty;                           %   quadtatic penalty on weights:  1/2 par.quadrpenpar*||vcwb||^2
	par.func_u=@(vcu,Z,par) err_nnfgh_u(vcu,Z,par,ytrain);    %   Discrepancy term:     1/2 ||Y- (v'Phi(U)+c)||^2
	par.multA= @(vcwb,par)  multWX(vcwb,par,Xtrain);          % user function   y=Ax
	par.multAt=@(vcu,par)   multUXt(vcu,par,Xtrain);          % user function  y=A'*x

if FlagSesop    %%%%% %     Optimize using SESOPTN tool       %%%%%%%%%%%%%
	
	options.report_func=@(x,report,par) nnreport(x,report,par,Xtrain,ytrain,Xtest,ytest);   % User function  to display  iteration progress;
	tic
	vcwb=sesoptn(vcwb_init,  par.func_u,  par.func_x, par.multA, par.multAt,options,par);  
	toc
end

if FlagFminunc  %%%%%       Optimize using matlab optimization toolbox function fminunc:   BFGS quasi-Newton
	par.flagXnew=1;
	%options=optimset('LargeScale','on','GradObj','on','Display','iter','MaxIter',1000);
	options=optimset('LargeScale','off','GradObj','on','Display','iter','MaxIter',200);
	tic
	vcwb=fminunc(  @(vcwb,par)err_nnfg(vcwb,par,Xtrain,ytrain),     vcwb_init,  options,par);  
	toc
end



% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% %
% %  Optimization using fminunc.m from Matlab Optimization toolbox 
% %
% %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%options1 =optimset('LargeScale','off','GradObj','on','Display','iter','MaxIter',200);
% disp('Optimization using fminunc.m from Matlab Optimization toolbox')
% c1 =fminunc( @fg, vcwb_init, options1, par);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%  Optimization using minFunc of Mark Schmidt 
%                              http://www.cs.ubc.ca/~schmidtm
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if FlagMinFunc
	disp('Optimization using L-BFGS: minFunc of Mark Schmidt')
	addpath('../../minFunc')
	par.flagXnew=1;
	options1 =optimset('GradObj','on','Display','iter','MaxIter',300);
	%options1.Display='final';
	options1.Method = 'lbfgs';
	options1.TolFun=1e-10;
	options1.TolX=1e-10;
	options1.Corr=8; % Number of previous steps used in L-BFGS update

	options1.MaxFunEvals=700;
	tic;
	[vcwb,f] =minFunc( @fg, vcwb_init, options1,par);
	toc
	fprintf('f_opt=%g ',f)
end

