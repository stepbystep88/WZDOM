function [x, fval, exitFlag, output] = bsPostInv1DTraceByDLSR(d, G, xInit, Lb, Ub, regParam, parampkgs, options)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code is designed for 1D seismic inversion using dictionary learning
% and sparse representation technique
%
% Programmed by: Bin She (Email: bin.stepbystep@gmail.com)
% Programming dates: Nov 2019
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
    if ~isfield(regParam, 'lambda') || regParam.lambda < 0
        % find the best regularization parameter
        regParam.lambda = bsFindBestRegParameter(options, inputObjFcnPkgs, xInit, Lb, Ub);
    end
           
    GBOptions = options.GBOptions;
    inputObjFcnPkgs{2, 3} = regParam.lambda;
    GBOptions.maxIter = options.innerIter;
    
    % create packages for sparse inversion 
    GSparseInvParam = bsInitDLSRPkgs(parampkgs, regParam.lambda, regParam.gamma, sampNum, G);
    
    ncell = GSparseInvParam.ncell;
    sizeAtom = GSparseInvParam.sizeAtom;
    
    midX = [];
    midF = [];
    
    for i = 1 : options.maxIter
        
        % change the current initial guess
        inputObjFcnPkgs{2, 2} = [];
        
%         if isfield(GSparseInvParam, 'isChangeK') && GSparseInvParam.isChangeK && i >= 3
%             GSparseInvParam.sparsity = 2;
%         end
        
        if isfield(GSparseInvParam, 'isScale') && GSparseInvParam.isScale
            Gx = mainData.A * xInit;
            c = (Gx' * mainData.B) / (Gx' * Gx);
            mainData.A = mainData.A * c;
            inputObjFcnPkgs{1, 2} = mainData;
        end
        
        [xOut, fval, exitFlag, output_] = bsGBSolveByOptions(inputObjFcnPkgs, xInit, Lb, Ub, GBOptions);
        
        if GBOptions.isSaveMiddleRes
            midX = [midX, output_.midResults.x];
            midF = [midF, output_.midResults.f];
        end
        
        % sparse reconstruction
        xI = exp(xOut);
        
        vps = zeros(sizeAtom, ncell);
        for j = 1 : ncell
            js = GSparseInvParam.index(j);
            vps(:, j) = xI(js : js+sizeAtom-1);
        end
        
        gammas = omp(GSparseInvParam.DIC'*vps, GSparseInvParam.omp_G, GSparseInvParam.sparsity);
        new_ips = GSparseInvParam.DIC *  gammas;
        
        %% reconstruct model by equations
        switch GSparseInvParam.reconstructType
            case 'equation'
                xNew = regParam.gamma * xI;
                % get reconstructed results by equation
                for j = 1 : ncell
                    xNew = xNew + GSparseInvParam.R{j}' * new_ips(:, j);
                end

                xNew = GSparseInvParam.invR * xNew;
            case 'simpleAvg'
                % get reconstructed results by averaging patches
                reconstruct_ips = bsAvgPatches(new_ips, GSparseInvParam.index, sampNum);
                xNew = reconstruct_ips * regParam.gamma + xI * (1 - regParam.gamma);
        end

        %% reconstruct model by 
        xNew = log(xNew);
        
        % projection
        xInit = bsPFunc(xNew, output_.options.bounds);
%         xInit = xNew;

%         figure;
%         plot(1:length(xNew), exp(xNew), 'r'); hold on;
%         plot(1:length(xNew), exp(xOut), 'g');
            
    end
    
    switch GSparseInvParam.isSparseRebuild
        case 1
            x = xInit;
        case 0 
            x = xOut;
        otherwise
            error('GSparseInvParam.isSparseRebuild must either 1 or 0. \n');
    end
    
    output.midResults.x = midX;
    output.midResults.f = midF;
    output.regParam = regParam;
    output.parampkgs = parampkgs;
    
end

function GSparseInvParam = bsInitDLSRPkgs(GSparseInvParam, lambda, gamma, sampNum, G)
    
    if isfield(GSparseInvParam, 'omp_G')
        return;
    end

    validatestring(string(GSparseInvParam.reconstructType), {'equation', 'simpleAvg'});
    validateattributes(gamma, {'double'}, {'>=', 0, '<=', 1});
    
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
%     GSparseInvParam.invG = inv(G' * G + gamma * tmp + lambda * eye(sampNum));
    GSparseInvParam.omp_G = GSparseInvParam.DIC' * GSparseInvParam.DIC;
    
end
