% SESOP-TN tool for Very Large scale Smooth Unconstrained Optimization
% implements the following methods:
% 
%     * Polak-Ribiere nonlinear Conjugate Gradients
%     * Truncated Newton
%     * SESOP:     Sequential Subspace Optimization
%     * SESOP-TN:  Sequential Subspace Optimization combined with Truncated Newton
%     * PCD,SSF,FISTA for L1-L2 optimization (see examples in Sesop_SPmagazine directory)
%
%
% Files:
%
% sesoptn.m            - function implementing SESOP-TN and other methods
% sesoptn_optionset.m  - setup options for sesoptn
%
%
% DIRECTORIES:
%
% Sesop_Basic_BP         - Basis Pursuit simple example
% Sesop_PCD_BasisPursuit - Basis Pursuit via PCD-SESOP (parallel coordinate descent) 
% Sesop_NeuralNet        - Neural Nnet Training. Applications: Signal Denoising and Prediction
% Sesop_SVM              - L1-L2 Support Vector Machine for Pattern Recognition 
% Sesop_SPmagazine       - Compressive Sensing, Image Deblurring and Tomography via L1-L2 optimization
% util_mz                - utilities: various service functions
% 
%
%
%
% References:
% 
% 1. Guy Narkiss and  Michael  Zibulevsky (2005). "Sequential Subspace Optimization Method for Large-Scale Unconstrained Problems", Tech. Report CCIT No 559, EE Dept., Technion. pdf file
% 2. Stephen G. Nash, “A Survey of Truncated-Newton Methods”, Journal of Computational and Applied Mathematics, 124 (2000), pp. 45-59.
% 3. Guy Narkiss and Michael Zibulevsky (2005). "Support Vector Machine via Sequential Subspace Optimization", Tech. Report CCIT No 557, EE Dept., Technion. pdf file
% 4. Michael Elad, Boaz Matalon, and Michael Zibulevsky, "Coordinate and Subspace Optimization Methods for Linear Least Squares with Non-Quadratic Regularization", Applied and Computational Harmonic Analysis, Vol. 23, pp. 346-367, November 2007. pdf file
% 5. Michael Zibulevsky and Michael Elad,"L1-L2 Optimization in Signal and Image Processing: Iterative Shrinkage and Beyond"
%    IEEE Signal Processing Magazine, to appear
%
% Michael Zibulevsky  06.08.2008; 05.02.2009; 26.10.2010
%
% Copyright (c) 2008-2010. All rights reserved. Free for academic use. No warranty 

