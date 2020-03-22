% L1-L2 Support Vector Machine for Pattern Recognition using SESOPTN - 
% Optimization tool for very large scale  smooth unconstrained minimization  
%
% Let matrix [X,y] contain training data points x_k and +-1 labels y_k at rows k=1,..,m
% We create matrix  A = [-diag(y)*X, -y];
%
% SVM training via smooth unconstrained optimization:
%
%            min_wb  mu*sum(abs_smooth(w)) + c*svm_penalty(A*[w;b])  
%
%
% See also:
% Narkiss, G. and Zibulevsky, M. (2005). "Support Vector Machine via Sequential 
% Subspace Optimization", Tech. Report CCIT No 557, EE Dept., Technion.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%   FILES:      
%
%   svm_sesoptn       - main script: L1-L2 Support Vector Machine 
%   svm_diag_penalty  - SVM penalty for separation vector w:           mu*sum(abs_smooth(w)) 
%   svm_penalty       - SVM penalty for violation of separating strip: mu1*sum(abs_smooth(z+1)+z+1)
%   svm_penalty_huber - SVM quadratic-linear penalty for violation of separating strip
%   svm_show_progress - Plot SVM error rate on validation set through SESOP ierations
%
%   uctptest_sesop    -  script for UCTP test optimization problems by Hans Bruun Nielsen
%   uctpval           -  function and gradient for UCTP problems
%   uctpval_sesop     -  SESOP interface for UCTP problems

% Michael Zibulevsky 10.07.2008;   05.08.2008; 03.02.2009
%
% Copyright (c) 2008 - 2009. All rights reserved. Free for academic use. No warranty 

