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
    GSParam = bsInitDLSRPkgs(parampkgs, regParam.gamma, sampNum);
    inline = options.inline;
    crossline = options.crossline;
    
    ncell = GSParam.ncell;
    sizeAtom = GSParam.sizeAtom;
    
    patches = zeros(sizeAtom*3+GSParam.nSpecialFeat, ncell);
    
    midX = [];
    midF = [];
    data = zeros(sampNum, 4);
    newData = data;
    maxIter = options.maxIter;
    lambda = regParam.lambda(1);
    
    
    for iter = 1 : maxIter
        
        % change the current initial guess
        inputObjFcnPkgs{2, 2} = [];
        if length(regParam.gamma) == 2
            gamma  = (maxIter - iter)*(regParam.gamma(2) - regParam.gamma(1))/(maxIter - 1) + regParam.gamma(1);
        else
            gamma = regParam.gamma;
        end
        
        if length(regParam.lambda) == 2
            lambda  = lambda * regParam.lambda(2);
            inputObjFcnPkgs{2, 3} = lambda;
        else
            inputObjFcnPkgs{2, 3} = regParam.lambda;
        end
        
        if length(options.initRegParam) == 2
            initLambda = initLambda * options.initRegParam(2);
            inputObjFcnPkgs{3, 3} = initLambda;
        end
        
        if iter == 1
            inputObjFcnPkgs{2, 3} = 0;
        end

        [xOut, fval, exitFlag, output_] = bsGBSolveByOptions(inputObjFcnPkgs, xInit, Lb, Ub, GBOptions);

        if GBOptions.isSaveMiddleRes
            midX = [midX, output_.midResults.x];
            midF = [midF, output_.midResults.f];
        end
        
        % sparse reconstruction
        [data(:, 2), data(:, 3), data(:, 4)] = bsPreRecoverElasticParam(xOut, mode, lsdCoef);
        
        for j = 1 : ncell
            js = GSParam.index(j);
            
            for i = 1 : 3
                
                sPos = sizeAtom*(i-1) + 1 + GSParam.nSpecialFeat;
                ePos = sPos + sizeAtom - 1;

                iData = data(js : js+sizeAtom-1, i+1);
                patches(sPos:ePos, j) = iData;
            end
        end
        
        % add location or/and time information
        if GSParam.trainDICParam.isAddLocInfo && GSParam.trainDICParam.isAddTimeInfo
            patches(1:GSParam.nSpecialFeat, :) = [ones(1, ncell) * inline; ones(1, ncell) * crossline; 1 : ncell];
        elseif GSParam.trainDICParam.isAddLocInfo
            patches(1:GSParam.nSpecialFeat, :) = [ones(1, ncell) * inline; ones(1, ncell) * crossline;];
        elseif GSParam.trainDICParam.isAddTimeInfo
            patches(1:GSParam.nSpecialFeat, :) = 1 : ncell;
        end
            
        % normalization
        normal_patches = (patches - GSParam.min_values) ./ (GSParam.max_values - GSParam.min_values);
        
        if strcmp(GSParam.trainDICParam.feature_reduction, 'all')
            normal_patches = GSParam.output.B' * normal_patches;
        end
    
        if GSParam.isModifiedDIC
            normal_patches = GSParam.M  * normal_patches;
        end
        
        gammas = omp(GSParam.DIC'*normal_patches, ...
                    GSParam.omp_G, ...
                    GSParam.sparsity);
        new_patches = GSParam.DIC *  gammas;
        
        if strcmp(GSParam.trainDICParam.feature_reduction, 'all')
            new_patches = GSParam.output.B * new_patches;
        end
        
        new_patches = new_patches .* (GSParam.max_values - GSParam.min_values) + GSParam.min_values;
        
        
        
        %% reconstruct model by equations
        for i = 1 : 3
            sPos = sizeAtom*(i-1) + 1 +  + GSParam.nSpecialFeat;
            ePos = sPos + sizeAtom - 1;
            
            i_new_patches = new_patches(sPos:ePos, :);
            switch GSParam.reconstructType
                case 'equation'
                    avgLog = gamma * data(:, i+1);
                    % get reconstructed results by equation
                    for j = 1 : ncell
                        
                        avgLog = avgLog + GSParam.R{j}' * i_new_patches(:, j);
                    end

                    newData(:, i+1) = GSParam.invR * avgLog;
                case 'simpleAvg'
                    % get reconstructed results by averaging patches
                    avgLog = bsAvgPatches(i_new_patches, GSParam.index, sampNum);
                    newData(:, i+1) = avgLog * gamma + data(:, i+1) * (1 - gamma);
            end
        end
        
        
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

function GSParam = bsInitDLSRPkgs(GSParam, gamma, sampNum)
    
    if isfield(GSParam, 'omp_G')
        return;
    end

    validatestring(string(GSParam.reconstructType), {'equation', 'simpleAvg'});
    validateattributes(gamma, {'double'}, {'>=', 0, '<=', 1});
    
    trainDICParam = GSParam.trainDICParam;
    
    sizeAtom = trainDICParam.sizeAtom;
    nAtom= trainDICParam.nAtom;
    GSParam.sizeAtom = sizeAtom;
    GSParam.nAtom = nAtom;
    GSParam.nrepeat = sizeAtom - GSParam.stride;
    
    GSParam.nSpecialFeat = trainDICParam.isAddLocInfo * 2 + trainDICParam.isAddTimeInfo;
    
    index = 1 : GSParam.stride : sampNum - sizeAtom + 1;
    if(index(end) ~= sampNum - sizeAtom + 1)
        index = [index, sampNum - sizeAtom + 1];
    end
    
    GSParam.index = index;
    GSParam.ncell = length(index);
    [GSParam.R] = bsCreateRMatrix(index, sizeAtom, sampNum);
   
    tmp = zeros(sampNum, sampNum);
    for iCell = 1 : GSParam.ncell
        tmp = tmp + GSParam.R{iCell}' * GSParam.R{iCell};
    end
    GSParam.invTmp = tmp;
    GSParam.invR = inv(gamma(1) * eye(sampNum) + GSParam.invTmp);
    
    if GSParam.isModifiedDIC
        I = eye(sizeAtom * 3);
        oneSa = ones(sizeAtom, sizeAtom);
        Z = zeros(sizeAtom, sizeAtom);
        cOne = {oneSa Z Z;
              Z oneSa Z;
              Z Z oneSa};
        GSParam.M = I + GSParam.a / sizeAtom * cell2mat(cOne);

        MDIC = GSParam.M * GSParam.DIC;
        % normalize the modified dictionary
        for j = 1 : size(MDIC, 2)
            MDIC(:, j) = MDIC(:, j) / norm(MDIC(:, j));
        end
        
        GSParam.DIC = MDIC;
        GSParam.omp_G = GSParam.MDIC' * GSParam.MDIC;
    else
        GSParam.omp_G = GSParam.DIC' * GSParam.DIC;
    end
        
    rangeCoef = GSParam.rangeCoef;
    
    GSParam.min_values = repmat(rangeCoef(:, 1), 1, GSParam.ncell);
	GSParam.max_values = repmat(rangeCoef(:, 2), 1, GSParam.ncell);
    
end


