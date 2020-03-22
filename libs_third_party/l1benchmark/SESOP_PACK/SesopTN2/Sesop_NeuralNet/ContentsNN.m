% Neural Nnet Training via SesopTN optimization: Signal denoising and prediction
%
% SESOPTN    Optimization tool for very large scale  smooth unconstrained minimization  
% 
%   FILES:      
%
%   sesoptn            - unconstrained optimization via Sequential subspace optimization (SESOP) 
%                        combined with Truncated Newton
%
%   sesoptn_optionset  - Setup options structure for sesoptn optimization function
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                    Training Neural Nnet via SesopTN optimization
%
%     min_vcWb  1/2 ||Ytrain-NN(Xtrain; vcWb)||^2 +  1/2 par.quadrpenpar||vcWb||^2
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% Scripts
%
%   nnSignalDenoise    - Training single layer Neural Net (NN) for signal  denoising using SESOPTN optimization tool
%   nnpredictSesop     - Training single layer Neural Net for time series prediction using SESOPTN optimization tool
%
% Functions
%   err_nnfg           - Error function of Single-layer feed-forward Neural Net and its gradient with respect to weights   
%   err_nnfgh_u        - Error function of  Neural Net with respevt to U=WX, its gradient and Hessian-vector (matrix) multiplication
%   multWX             - vcub=multA(vcwb,par):  U=W*X; vcub=[v;c;U(:);b];
%   multUXt            - Adjoint operator to multWX(): W=U*X'; vcwb=[v;c;W(:);b];
%   nnet               - Simgle layer Neural Net:  y=v'*sigmoid(W*X+b*1')+c
%   nnreport           - Plot Neural Net error during learning SESOP_TN iterations
%   sigmoid_mz         - Sigmoid function:   f(t) = t/(1+|t|);  f_eps(t)=f(t/eps);

% Michael Zibulevsky,  04.08.2008 03.02.2009
%
% Copyright (c) 2008-2009. All rights reserved. Free for academic use. No warranty 
