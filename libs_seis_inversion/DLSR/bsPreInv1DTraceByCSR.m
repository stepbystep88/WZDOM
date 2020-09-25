function [x, fval, exitFlag, output] = bsPreInv1DTraceByCSR(d, G, xInit, Lb, Ub, regParam, parampkgs, options, mode, lsdCoef)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code is designed for prestack 1D seismic inversion using dictionary
% learning and collaboration sparse representation representation
%
% Programmed by: Bin She (Email: bin.stepbystep@gmail.com)
% Programming dates: Dec 2019
% -------------------------------------------------------------------------
% Input
% 
% d: observed seismic data, a vector.
% 
% G: forward operator made up by the wavelet information
% 
% xInit: initial guess of model parameters
% 
% Lb: lower boundary of x
% 
% Ub: upper boundary of x
% 
% regParam: regularization parameter. If it is empty, I will start a search
% process to find the optimal regParam.
%
% parampkgs: specail parameter for diffrent methods.
% 
% options: options parameters for 1D seismic inversion. See function
% bsCreateSeisInv1DOptions
% 
% regFunc: regularization function handle. It could be @bsReg*.m
%
% -------------------------------------------------------------------------
% Output
%
% x          is a column vector; refers to the estimated result.
%
% fval          the objective function value at the last iteration
%
% exitFlag      corresponds to the stopCriteria
% see function bsCheckStopCriteria
%
% output    a struct, in general, it has
% output.iterations: the number of iterations
% output.nfev: the number of function evaluations
% ouput.midResults: the middle results during the iteration process
% output.regParam: the regularization parameters used
% -------------------------------------------------------------------------

    % create mainData
    mainData.A = G;
    mainData.B = d;
    sampNum = length(xInit)/3;
    initLambda = options.initRegParam(1);
    
    % re-organize the input objective function pakages
    inputObjFcnPkgs = {
        options.mainFunc,       mainData,   1; 
        @bsReg1DTKInitModel,    struct('xInit', xInit), initLambda;
        @bsReg1DTKInitModel,    struct('xInit', xInit), initLambda;
    };
    
    % if the regParam is not given, I search it by a search subroutine
    % which is save in options.searchRegParamFcn. 
    if ~isfield(regParam, 'lambda')
        % find the best regularization parameter
        regParam.lambda = bsFindBestRegParameter(options, inputObjFcnPkgs, xInit, Lb, Ub);
    end
           
    GBOptions = options.GBOptions;
    inputObjFcnPkgs{2, 3} = regParam.lambda;
    GBOptions.maxIter = options.innerIter;
    
    % create packages for sparse inversion 
%     GSParam = bsInitDLSRPkgs(parampkgs, regParam.gamma, sampNum);
    GSParam = bsInitGSparseParam(parampkgs, sampNum, 1);
    inline = options.inline;
    crossline = options.crossline;
    gamma = regParam.gamma(1);
    
%     ncell = GSParam.ncell;
%     sizeAtom = GSParam.sizeAtom;
    
%     patches = zeros(sizeAtom*3+GSParam.nSpecialFeat, ncell);
    
    midX = [];
    midF = [];
    data = zeros(sampNum, 4);
%     newData = data;
    maxIter = options.maxIter;
%     lambda = regParam.lambda(1);
    
    nDic = (size(GSParam.DIC, 1) - GSParam.nSpecialFeat) / GSParam.sizeAtom;
    
    for iter = 1 : maxIter
        
        % change the current initial guess
        inputObjFcnPkgs{2, 2} = [];
        inputObjFcnPkgs{2, 3} = regParam.lambda;

        [xOut, fval, exitFlag, output_] = bsGBSolveByOptions(inputObjFcnPkgs, xInit, Lb, Ub, GBOptions);

        if GBOptions.isSaveMiddleRes
            midX = [midX, output_.midResults.x];
            midF = [midF, output_.midResults.f];
        end
        
        % sparse reconstruction
        [data(:, 2), data(:, 3), data(:, 4)] = bsPreRecoverElasticParam(xOut, mode, lsdCoef);
        
        % 接下来的代码时对模型m进行稀疏重构
        vp_vs = data(:, 2)./data(:,3);
        
        switch nDic
            case 3
                [out, gammas] = bsSparseRebuildOneTrace(GSParam, {data(:, 2), data(:, 3), data(:, 4)}, gamma, inline, crossline);
            case 4
                [out, gammas] = bsSparseRebuildOneTrace(GSParam, {data(:, 2), data(:, 3), data(:, 4), vp_vs}, gamma, inline, crossline);
        end
            
        
        
        newData = [data(:, 1), out(:, 1:3)];
        
        %% reconstruct model by 
        xInit = bsPreBuildModelParam(newData, mode, lsdCoef);
        
    end
    
    switch GSParam.isSparseRebuild
        case 1
            x = xInit;
        case 0 
            x = xOut;
        otherwise
            error('GSParam.isSparseRebuild must either 1 or 0. \n');
    end
    
    output.midResults.x = midX;
    output.midResults.f = midF;
    output.regParam = regParam;
    output.parampkgs = GSParam;
    
    output.gammas = gammas;
end


