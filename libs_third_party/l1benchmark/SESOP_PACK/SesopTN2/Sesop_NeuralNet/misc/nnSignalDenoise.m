% Training single layer Neural Net (NN) for signal  denoising using SESOPTN optimization tool

% Michael Zibulevsky, 24.07.2008; 25.07.08; 27.07.08; 31.07.08; 01.08.08; 03.02.2009
%
% Copyright (c) 2008 - 2009. All rights reserved. Free for academic use. No warranty 

%profile off
%profile on
%addpath('../');addpath('../util_mz');
addpath('../');addpath('../..');addpath('../../util_mz');

global GlobalErrTrain GlobalErrTest GlobalNNiter
GlobalErrTrain=[]; GlobalErrTest =[]; GlobalNNiter=[];


INIT_SIGNALS=0;          % 1 - generate new signals; 0 - use signals from the previous run  
INIT_DATA_MATRIX=1;      % 1 - generate new data matrix ;  0 - use data matrix from the previous run  
INIT_VARIABLES=0;        % 1 - generate random starting point; 0 - use old starting point
continue_old  =0;        % 1- run again starting from last achieved solution with possibly changed  parameters

%%%%% Optimize with: 
FlagSesop =   1;   %  1 - sesoptn tool;   
FlagFminunc = 0;   %  1 - fminunc BFGS Quasi-Newton (0ptimization toolbox of MATAB)
FlagMinFunc = 0;   %  1 - minFunc L-BFGS Quasi-Newton


options=sesoptn_optionset;  % Get default options structure (see comments in sesoptn_optionset.m) 
par.eps_sigmoid=0.7;        % Sigmoid parameter for sigmoid_mz function
par.quadrpenpar=1e-1;       % Quadratic penalty for NN weights: used for regularization, to avoid overfitting
par.nneurons=5;             % Number of neurons

par.nnreport_figure_handle= figure('Position',[10 40 400 600], 'Name',sprintf('NN error function: penpar=%g',par.quadrpenpar));



T=3000;           % Signal length:   the same for training and test signal
delta=20;         % Half-interval around the reconstructed sample to be feed  into NN;
sigma_noise=0.5;  % Noise std in simulated data


% options.TESTGRAD=1;  % 1 - Test gradient numerically, plot testing results (without run optimization); 0 - don't test
%par.nneurons=2; delta=1; T=10; %Small problem for gradient testing


if INIT_SIGNALS,      INIT_DATA_MATRIX=1; INIT_VARIABLES=1;   continue_old  =0;    end
if INIT_DATA_MATRIX,  INIT_VARIABLES=1;continue_old  =0;                           end
if continue_old,      INIT_SIGNALS=0;  INIT_DATA_MATRIX=0; INIT_VARIABLES=0;       end



	
	
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%      Generate clean and noisy signals (bumps + sinusoid)  for training and test sets
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	

if INIT_SIGNALS,  
	
	tt=[1:T]'/T;  % Time running from 0 to 1
	
	S00Train= cumsum(full(sprandn(T,1,0.02)))  +  3* sin(20*tt); % bumps + sinusoid
	NoiseTrain=sigma_noise*randn(T,1);
	SnoisyTrain=S00Train+NoiseTrain;
	
	S00Test=  cumsum(full(sprandn(T,1,0.02)))  +  2* sin(10*tt).^2;
	NoiseTest=sigma_noise*randn(T,1);
	SnoisyTest=S00Test+NoiseTest;
end

	

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                   Creating  training and test set data matrices:
%
%                          Xtrain, ytrain;   Xtest, ytest
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



if INIT_DATA_MATRIX,

	disp('Creating  training and test data matrices');

    K=T-2*delta;           % Number of training examples
	L= 2*delta+1;           % Interval around the reconstructed sample to be feed  into NN;
		
	%%%%%%%%     Training set    %%%%%

	Xtrain=zeros(L,K);
	ytrain=zeros(1,K);
	
	k=0;
	for t=1+delta: T-delta,   
		k=k+1;
		Xtrain(1:L,k)=SnoisyTrain(t-delta : t+delta) - SnoisyTrain(t);
		ytrain(k)=NoiseTrain(t);
	end

	%%%%%%%     Test set   %%%%%
	
	Xtest=zeros(L,K);
	ytest=zeros(1,K);

	k=0;
	for t=1+delta: T-delta,  
		k=k+1;
		Xtest(1:L,k)=SnoisyTest(t-delta : t+delta) - SnoisyTest(t);
		ytest(k)=NoiseTest(t);
	end

end 




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%           Random initialization of Neutal Net (NN)  weights  v,W,c,b  
%
%                where     NN(X; vcwb)  =  v'phi(WX+b*1') + c 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if INIT_VARIABLES   

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

if continue_old,	vcwb_init=vcwb;      %  Skip weights initialization, if you want just improve previous optimization results
else                vcwb_init=vcwb0;
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                Learning of Neural Nnet  weights via optimization:   
%
%     min_vcwb  1/2 ||Ytrain - NN(Xtrain;vcwb)||^2 +  1/2 par.quadrpenpar*||vcwb||^2
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




fprintf('\n Perform learning of Neural Net  weights via optimization:\n \n ');

	par.func_x=@diag_quadr_penalty;                                    %   quadtatic penalty on weights:  1/2 par.quadrpenpar*||vcwb||^2
	par.func_u=@(vcu,Z,par) err_nnfgh_u(vcu,Z,par,ytrain);    %   Discrepancy term:     1/2 ||Y- (v'Phi(U)+c)||^2
	par.multA= @(vcwb,par)  multWX(vcwb,par,Xtrain);             % user function   y=Ax
	par.multAt=@(vcu,par)   multUXt(vcu,par,Xtrain);              % user function  y=A'*x

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
	options1 =optimset('GradObj','on','Display','iter','MaxIter',1000);
	%options1.Display='final';
	options1.Method = 'lbfgs';
	options1.TolFun=1e-10;
	options1.TolX=1e-10;
	options1.Corr=8; % Number of previous steps used in L-BFGS update

	options1.MaxFunEvals=700;
	tic;
	[wb,f] =minFunc( @fg, vcwb_init, options1,par);
	toc
	fprintf('f_opt=%g ',f)
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%           Compute denoised signals  and plot results
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      


ind=1+delta:T-delta;
SdenoisedTrain = SnoisyTrain(ind) - nnet(vcwb,par,Xtrain)';  
SdenoisedTest  = SnoisyTest(ind)  - nnet(vcwb,par,Xtest)';

figure; subplot(311);
plot( [S00Train(ind); S00Test(ind)]  );grid; 
title('Original  signals:  Left half - Training, Right half - Test')
subplot(312);
plot([ [S00Train(ind); S00Test(ind)]  [SnoisyTrain(ind); SnoisyTest(ind)] ] );grid; 
title('Noisy signals')
subplot(313);
plot([ [S00Train(ind); S00Test(ind)]  [SdenoisedTrain;SdenoisedTest ]] );grid; 
title('Green - Recovered with Neural Network; Blue - Original')


%profile viewer