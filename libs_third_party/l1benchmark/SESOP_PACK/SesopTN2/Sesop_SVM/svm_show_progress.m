function show_progressSVM(w,report,par),
%Plot SVM error rate on validation set 

% Michael Zibulevsky, 05.08.2008        
%
% Copyright (c) 2008. All rights reserved. Free for academic use. No warranty 

global Global_nn_errors  Global_sesop_iter

Global_sesop_iter=[Global_sesop_iter report.Niter];

b=w(end);
w=w(1:end-1);
n=length(w);
m=length(par.y_validation);
%error_vector=( w'*par.X_validation + b*ones(1,m) );
error_vector=( w'*par.X_validation + b*ones(1,m) ).*par.y_validation' <0;
n_errors= sum( error_vector(:));
fprintf('num_errors %d \n',n_errors);
Global_nn_errors=[Global_nn_errors n_errors];
figure(par.report_figure_handle);subplot(211);plot(Global_sesop_iter,Global_nn_errors/m);
title('Frequency of errors in validation set');grid
subplot(212);plot(w);title('w - separating vector');grid