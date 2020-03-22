%Support Vector Machine using sesopmz

  
%global GlobalNgrads GlobalGradNorms  GlobalNfuncs GlobalFuncValues GlobalNmultA

%global DataMatrix
global nn_errors

INIT_DATA=1;  

options=sesoptn_optionset;  % Get default options structure (see comments in optionset_sesoptn.m) 

options.sesop_figure_name=sprintf('SESOPtn progress,   %d CG steps per TN iteration',options.max_iter_CGinTN);
options.report_func=@svm_show_progress;
par.report_figure_handle= figure('Position',[10 40 400 600], 'Name',sprintf('SVM validation error rate'));


n=200; % data dimension 
m=500; % number of training examples

FlagSubsets=0;
par.Nsubsets=1e20;

par.eps_smooth_abs=0.2;  % Smoothing parameter of absolute value approximation
par.c=3*n;               % Weight of penalty for violation of the separating band




% Generate data for classification:  Two Gaussian clouds shifted by delta in each coodinate
if INIT_DATA,
	delta=4/sqrt(n);
	m1=m/4;    X1=randn(n,m1); y1=ones(m1,1);
	m2=m-m1; X2=randn(n,m2) + delta*ones(n,m2); y2=-ones(m2,1);
	X=[X1 X2]; y=[y1;y2];
	
	ind=randperm(m);
	X=X(:,ind);y=y(ind);
	
	DataMatrix= [ -mulmd(X,y)', -y];
	w0=randn(n+1,1);

	X1=randn(n,m/2); y1=ones(m/2,1);
	X2=randn(n,m/2) + delta*ones(n,m/2); y2=-ones(m/2,1);
	par.X_validation=[X1 X2]; par.y_validation=[y1;y2];
	
end

%%%%%%%%%%%%  For debugging
%ind1=[5: m];                                   
%DataMatrix= [mulmd(X,y)', y];
%DataMatrix= [-mulmd(X(:,ind1),y(ind1))', -y(ind1)];
%DataMatrix=DataMatrix(ind1,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%figure;imagesc(x00);colorbar;return





nn_errors = [];

par.func_u = @svm_penalty;
par.func_x = @svm_diag_penalty;

if FlagSubsets
	%par.figreport_num=989;
	par.multA= @mult_global_matrix_subset;
	par.multAt=@mult_global_matrix_subset_adj;
	[w,report]=sesop_subset(w0(:),@svm_penalty,@svm_diag_penalty, par.multA, par.multAt, options,par);
else
	%par.figreport_num=987;
	par.multA= @(x,par)  multMatr(DataMatrix,x);     % user function   y=Ax
	par.multAt=@(x,par)  multMatrAdj(DataMatrix,x);  % user function  y=A'*x

	%par.multA= @mult_global_matrix;
	%par.multAt=@mult_global_matrix_adj;

	%par.multA= @mult_global_matrix_subset;
	%par.multAt=@mult_global_matrix_subset_adj;
	
	%[w,report]=sesopmz(w0(:),@svm_penalty,@svm_diag_penalty, par.multA, par.multAt, options,par);
	[w,report]=sesoptn(w0(:),par.func_u, par.func_x, par.multA, par.multAt,options,par);

end

b=w(end);
w=w(1:end-1);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            Show results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

