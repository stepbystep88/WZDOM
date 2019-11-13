function [x, fval, exitFlag, output] = bsPostInv1DTraceByDLSR(d, G, xInit, Lb, Ub, regParam, parampkgs, options)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code is designed for 1D seismic inversion using regularization
% technique
%
% Programmed by: Bin She (Email: bin.stepbystep@gmail.com)
% Programming dates: May 2019
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
    sampNum = length(xInit);
    
    % re-organize the input objective function pakages
    inputObjFcnPkgs = {
        options.mainFunc,       mainData,   1; 
        @bsReg1DTKInitModel,    struct('xInit', xInit), options.initRegParam;
        @bsReg1DTKInitModel,    struct('xInit', xInit), options.initRegParam;
    };
    
    % if the regParam is not given, I search it by a search subroutine
    % which is save in options.searchRegParamFcn. 
    if ~isfield(regParam, 'lambda') || regParam.lambda <= 0
        % find the best regularization parameter
        regParam.lambda = bsFindBestRegParameter(options, inputObjFcnPkgs, xInit, Lb, Ub);
    end
           
    GBOptions = options.GBOptions;
    inputObjFcnPkgs{2, 3} = regParam.lambda;
    GBOptions.maxIter = options.innerIter;
    
    % create packages for sparse inversion 
    parampkgs = bsInitDLSRPkgs(parampkgs, regParam.lambda, regParam.gamma, sampNum, G);
    GSparseInvParam = parampkgs.GSparseInvParam;
    
    ncell = GSparseInvParam.ncell;
    sizeAtom = GSparseInvParam.sizeAtom;
    
    midX = [];
    midF = [];
    
    for i = 1 : options.maxIter
        
        % change the current initial guess
        inputObjFcnPkgs{2, 2} = [];
        [xOut, fval, exitFlag, output_] = bsGBSolveByOptions(inputObjFcnPkgs, xInit, Lb, Ub, GBOptions);
        
        if GBOptions.isSaveMiddleRes
            midX = [midX, output_.midResults.x];
            midF = [midF, output_.midResults.f];
        end
        
        % sparse reconstruction
        xI = exp(xOut);
        
        xNew = regParam.gamma * xI;
        vps = zeros(sizeAtom, ncell);
        for j = 1 : ncell
            js = GSparseInvParam.index(j);
            vps(:, j) = xI(js : js+sizeAtom-1);
        end
        
        gammas = omp(GSparseInvParam.DIC'*vps, GSparseInvParam.omp_G, GSparseInvParam.sparsity);
        new_vps = GSparseInvParam.DIC *  gammas;
        
        for j = 1 : ncell
            xNew = xNew + GSparseInvParam.R{j}' * new_vps(:, j);
        end
        
        xNew = GSparseInvParam.invR * xNew;
        xNew = log(xNew);
        
        % projection
        xInit = bsPFunc(xNew, output_.options.bounds);
%         xInit = xNew;

            
    end
    
    x = xInit;
    
    output.midResults.x = midX;
    output.midResults.f = midF;
    output.regParam = regParam;
    output.parampkgs = parampkgs;
    
end

function parampkgs = bsInitDLSRPkgs(parampkgs, lambda, gamma, sampNum, G)
    
    if isfield(parampkgs.GSparseInvParam, 'invR')
        return;
    else
        GSparseInvParam = parampkgs.GSparseInvParam;
    end

    
    [sizeAtom, nAtom] = size(GSparseInvParam.DIC);
    GSparseInvParam.sizeAtom = sizeAtom;
    GSparseInvParam.nAtom = nAtom;
    GSparseInvParam.nrepeat = sizeAtom - GSparseInvParam.stride;
    
    index = 1 : GSparseInvParam.stride : sampNum - sizeAtom + 1;
    if(index(end) ~= sampNum - sizeAtom + 1)
        index = [index, sampNum - sizeAtom + 1];
    end
    
    GSparseInvParam.index = index;
    GSparseInvParam.ncell = length(index);
    [GSparseInvParam.R] = bsCreateRMatrix(index, sizeAtom, sampNum);
   
    tmp = zeros(sampNum, sampNum);
    for iCell = 1 : GSparseInvParam.ncell
        tmp = tmp + GSparseInvParam.R{iCell}' * GSparseInvParam.R{iCell};
    end
    GSparseInvParam.invTmp = tmp;
    GSparseInvParam.invR = inv(gamma * eye(sampNum) + GSparseInvParam.invTmp);
    GSparseInvParam.invG = inv(G' * G + gamma * tmp + lambda * eye(sampNum));
    GSparseInvParam.omp_G = GSparseInvParam.DIC' * GSparseInvParam.DIC;
    
    parampkgs.GSparseInvParam = GSparseInvParam;
end

