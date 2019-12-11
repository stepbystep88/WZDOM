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
    parampkgs = bsInitDLSRPkgs(parampkgs, regParam.gamma, sampNum);
    GSparseInvParam = parampkgs.GSparseInvParam;
    
    ncell = GSparseInvParam.ncell;
    sizeAtom = GSparseInvParam.sizeAtom;
    rangeCoef = GSparseInvParam.rangeCoef;
    patches = zeros(sizeAtom*3, ncell);
    
    midX = [];
    midF = [];
    data = zeros(sampNum, 4);
    newData = data;
    
    
    
    for iter = 1 : options.maxIter
        
        % change the current initial guess
        inputObjFcnPkgs{2, 2} = [];
        [xOut, fval, exitFlag, output_] = bsGBSolveByOptions(inputObjFcnPkgs, xInit, Lb, Ub, GBOptions);
        
        if GBOptions.isSaveMiddleRes
            midX = [midX, output_.midResults.x];
            midF = [midF, output_.midResults.f];
        end
        
        % sparse reconstruction
        [data(:, 2), data(:, 3), data(:, 4)] = bsPreRecoverElasticParam(xOut, mode, lsdCoef);
        
        for j = 1 : ncell
            js = GSparseInvParam.index(j);
            for i = 1 : 3
                sPos = sizeAtom*(i-1) + 1;
                ePos = sPos + sizeAtom - 1;
                
                iData = data(js : js+sizeAtom-1, i+1);
                % normalization
                patches(sPos:ePos, j) = (iData - rangeCoef(i, 1))/(rangeCoef(i, 2) - rangeCoef(i, 1));
            end
        end

        if GSparseInvParam.isModifiedDIC
            patches = GSparseInvParam.M  * patches;
        end
        
        gammas = omp(GSparseInvParam.DIC'*patches, ...
                    GSparseInvParam.omp_G, ...
                    GSparseInvParam.sparsity);
        new_patches = GSparseInvParam.DIC *  gammas;
        
        %% reconstruct model by equations
        for i = 1 : 3
            sPos = sizeAtom*(i-1) + 1;
            ePos = sPos + sizeAtom - 1;
            i_new_patches = new_patches(sPos:ePos, :) * (rangeCoef(i, 2) - rangeCoef(i, 1)) + rangeCoef(i, 1);
            
            switch GSparseInvParam.reconstructType
                case 'equation'
                    avgLog = regParam.gamma * data(:, i+1);
                    % get reconstructed results by equation
                    for j = 1 : ncell
                        
                        avgLog = avgLog + GSparseInvParam.R{j}' * i_new_patches(:, j);
                    end

                    newData(:, i+1) = GSparseInvParam.invR * avgLog;
                case 'simpleAvg'
                    % get reconstructed results by averaging patches
                    avgLog = bsAvgPatches(i_new_patches, GSparseInvParam.index, sampNum);
                    newData(:, i+1) = avgLog * regParam.gamma + data(:, i+1) * (1 - regParam.gamma);
            end
        end
        
        
        %% reconstruct model by 
        xInit = bsPreBuildModelParam(newData, mode, lsdCoef);
        
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

function parampkgs = bsInitDLSRPkgs(parampkgs, gamma, sampNum)
    
    if isfield(parampkgs.GSparseInvParam, 'omp_G')
        return;
    else
        GSparseInvParam = parampkgs.GSparseInvParam;
    end

    validatestring(string(GSparseInvParam.reconstructType), {'equation', 'simpleAvg'});
    validateattributes(gamma, {'double'}, {'>=', 0, '<=', 1});
    
    [sizeAtom, nAtom] = size(GSparseInvParam.DIC);
    sizeAtom = sizeAtom / 3;
    
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
    
    if GSparseInvParam.isModifiedDIC
        I = eye(sizeAtom * 3);
        oneSa = ones(sizeAtom, sizeAtom);
        Z = zeros(sizeAtom, sizeAtom);
        cOne = {oneSa Z Z;
              Z oneSa Z;
              Z Z oneSa};
        GSparseInvParam.M = I + GSparseInvParam.a / sizeAtom * cell2mat(cOne);

        MDIC = GSparseInvParam.M * GSparseInvParam.DIC;
        % normalize the modified dictionary
        for j = 1 : size(MDIC, 2)
            MDIC(:, j) = MDIC(:, j) / norm(MDIC(:, j));
        end
        
        GSparseInvParam.DIC = MDIC;
        GSparseInvParam.omp_G = GSparseInvParam.MDIC' * GSparseInvParam.MDIC;
    else
        GSparseInvParam.omp_G = GSparseInvParam.DIC' * GSparseInvParam.DIC;
    end
        
    parampkgs.GSparseInvParam = GSparseInvParam;
end


