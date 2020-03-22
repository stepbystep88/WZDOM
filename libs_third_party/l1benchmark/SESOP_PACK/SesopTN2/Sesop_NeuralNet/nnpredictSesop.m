% Training single layer Neural Net for time series prediction using SESOPTN optimization tool

% Michael Zibulevsky, 24.07.2008; 25.07.2008;   04.08.2008 
%
% Copyright (c) 2008. All rights reserved. Free for academic use. No warranty 

%profile off
%profile on
addpath('../');addpath('../util_mz');

INIT_DATA=1;        % 1 - generate data 
INIT_VARIABLES=0;   % 1 - generate random starting point
continue_old  =0;   % 1- run again starting from  achieved solution with possibly changed  parameters

if INIT_DATA,     INIT_VARIABLES=1;   continue_old=0;    end
if continue_old,  INIT_VARIABLES=0;   INIT_DATA=0;         end

%%%%% Optimize with: 
FlagSesop =1;    %  1 - sesoptn tool;   
FlagFminunc=0;   %  1 - fminunc BFGS Quasi-Newton (0ptimization toolbox of MATAB)
FlagMinFunc=0;   %  1 - minFunc L-BFGS Quasi-Newton

options=sesoptn_optionset;  % Get default options structure (see comments in optionset_sesoptn.m) 
options.sesop_figure_name=sprintf('SESOPtn progress,   %d CG steps per TN iteration',options.max_iter_CGinTN);


par.eps_sigmoid=0.7;
par.quadrpenpar=2e-2;
par.nneurons=3;

par.nnreport_figure_handle = figure('Position',[10 40 400 600], 'Name',sprintf('NN error function: penpar=%g',par.quadrpenpar));




Tmem=20; % Memory of the predictor;
%Tforc=500; % Forecast horizon
Tforc=1; % Forecast horizon

T=1024; % Train + Test overall time interval

%par.nneurons=2; Tmem=3; T=40;  % Small problem  for debugging
%options.TESTGRAD=1;  % 1 - Test gradient numerically,


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%      Generate signals to be predicted  for training and test sets
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


tt=[1:T];
if INIT_DATA,  
	%tmp=full(sprandn(T,1,0.02)); s=real(fft(tmp))';
	s=sin(0.5*tt) + cos(0.12*tt)+ cos(0.05*tt)+sin(0.055*tt); % +0.2*randn(size(tt));
	%s=load('dowdata.txt'); s=s(:,2)';s=diff(s);s=s/max(abs(s));
	% s=load('ta25data.txt'); s=flipud(s(:,4))'; s=diff(s);s=s/max(abs(s));
	
	%figure;plot(s); grid; return
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                   Creating  training and test set data matrices:
%
%                          Xtrain, ytrain;   Xtest, ytest
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



if INIT_VARIABLES
	T=length(s);
	Ttrain=floor(0.5*T);
	Ttest= floor(0.5*T);

	Xtrain=zeros(Tmem, Ttrain-Tmem-Tforc+1);
	Xtest =zeros(Tmem, Ttest-Tmem-Tforc+1);

	for i=1:Tmem
		Xtrain(i,:)=s(i:Ttrain-Tmem-Tforc+i);
	end
	ytrain=s(Tmem+Tforc : Ttrain);

	for i=1:Tmem
		Xtest(i,:)=s(Ttrain+i : Ttest+Ttrain-Tmem-Tforc+i);
	end
	ytest=s(Ttrain+Tmem+Tforc : Ttrain+Ttest);
	%figure;plot([ytrain ytest])

	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%           Random initialization of Neutal Net (NN)  weights  v,W,c,b  
%
%                where     NN(X; vcwb)  =  v'phi(WX+b*1') + c 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	
	[N,K]=size(Xtrain);
	par.Ktrain_samples=K;

	M=par.nneurons;
	
	v0=(1/sqrt(K))*rand(M,1);
	c0=0;
	W0=(1/sqrt(N))*randn(M,N);
	b0=0.1*randn(M,1);
	vcwb0=[v0;c0;W0(:);b0];
	vcwb=vcwb0;

end  


if continue_old,   vcwb_init=vcwb;
else               vcwb_init=vcwb0;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                Learning of Neural Nnet  weights via optimization:   
%
%     min_vcw  1/2 ||Ytrain -  NN(Xtrain; vcw)  ||^2 +  1/2 par.quadrpenpar*||vcw||^2
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n Perform learning of Neural Net  weights via optimization:\n \n ');

par.func_x=@diag_quadr_penalty;  %   quadtatic penalty on weights:  1/2 par.quadrpenpar*||vcwb||^2
par.func_u=@(vcu,Z,par) err_nnfgh_u(vcu,Z,par,ytrain);  %   Discrepancy term:     1/2 ||Y- (v'Phi(U)+c)||^2
par.multA= @(vcwb,par)  multWX(vcwb,par,Xtrain);        % user function   y=Ax
par.multAt=@(vcu,par)   multUXt(vcu,par,Xtrain);        % user function  y=A'*x

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
%           Compute predicted signals  and plot results
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 

ynn=nnet(vcwb,par,Xtrain);
figure;plot([ytrain' ynn' ]);grid; %plot([ytrain' ynn' sign(ytrain.*ynn)' ]);
title(sprintf('Training set: success ratio %g', mean( 0.5*sign(ytrain.*ynn)' + 0.5)))

ynn=nnet(vcwb,par,Xtest);
figure; plot([ytest' ynn' ]);grid;title(sprintf('Test: success ratio %g', mean( 0.5*sign(ytest.*ynn)' + 0.5)));

%profile viewer