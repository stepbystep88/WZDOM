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
    
    inline = options.inline;
    crossline = options.crossline;
    
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
    GSParam = bsInitDLSRPkgs(parampkgs, regParam.lambda, regParam.gamma, sampNum, G);
    
    ncell = GSParam.ncell;
    sizeAtom = GSParam.sizeAtom;
    
    midX = [];
    midF = [];
    vps = zeros(sizeAtom+GSParam.nSpecialFeat, ncell);
    
    for i = 1 : options.maxIter
        
        % change the current initial guess
        inputObjFcnPkgs{2, 2} = [];
        
%         if isfield(GSParam, 'isChangeK') && GSParam.isChangeK && i >= 3
%             GSParam.sparsity = 2;
%         end
%         if i == 1
%             inputObjFcnPkgs{2, 3} = 0;
%         else
%             inputObjFcnPkgs{2, 3} = regParam.lambda;
%         end
        
        if isfield(GSParam, 'isScale') && GSParam.isScale
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
        
        
        for j = 1 : ncell
            js = GSParam.index(j);
            vps(GSParam.nSpecialFeat+1 : end, j) = xI(js : js+sizeAtom-1);
        end
        
        % add location or/and time information
        if GSParam.trainDICParam.isAddLocInfo && GSParam.trainDICParam.isAddTimeInfo
            vps(1:GSParam.nSpecialFeat, :) = [ones(1, ncell) * inline; ones(1, ncell) * crossline; 1 : ncell];
        elseif GSParam.trainDICParam.isAddLocInfo
            vps(1:GSParam.nSpecialFeat, :) = [ones(1, ncell) * inline; ones(1, ncell) * crossline;];
        elseif GSParam.trainDICParam.isAddTimeInfo
            vps(1:GSParam.nSpecialFeat, :) = 1 : ncell;
        end
        
        switch GSParam.trainDICParam.normalizationMode
            case {'off', 'none'}
                
            otherwise
                % normalization
                vps = (vps - GSParam.min_values) ./ (GSParam.max_values - GSParam.min_values);
        end
        
        if ~strcmp(GSParam.trainDICParam.feature_reduction, 'off')
            reduced_vps = GSParam.output.B' * vps;
        else
            reduced_vps = vps;
        end
        
        gammas = omp(GSParam.DIC'*reduced_vps, GSParam.omp_G, GSParam.sparsity);
%         [gammas, oldGammas] = bsOMP(GSParam.DIC, vps, GSParam.omp_G, GSParam.sparsity, GSParam.neiborIndecies);
        new_ips = GSParam.DIC *  gammas;
        
        if strcmp(GSParam.trainDICParam.feature_reduction, 'all')
            new_ips = GSParam.output.B * new_ips;
        end
        
        switch GSParam.trainDICParam.normalizationMode
            case {'off', 'none'}
            otherwise
                % normalization
                new_ips = new_ips .* (GSParam.max_values - GSParam.min_values) + GSParam.min_values;
        end
        
        %% reconstruct model by equations
        switch GSParam.reconstructType
            case 'equation'
                xNew = regParam.gamma * xI;
                % get reconstructed results by equation
                for j = 1 : ncell
                    xNew = xNew + GSParam.R{j}' * new_ips(:, j);
                end

                xNew = GSParam.invR * xNew;
            case 'simpleAvg'
                % get reconstructed results by averaging patches
                reconstruct_ips = bsAvgPatches(new_ips(GSParam.nSpecialFeat+1:end,:), GSParam.index, sampNum);
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
    
end

function GSParam = bsInitDLSRPkgs(GSParam, lambda, gamma, sampNum, G)
    
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
    
    index = 1 : GSParam.stride : sampNum - sizeAtom + 1;
    if(index(end) ~= sampNum - sizeAtom + 1)
        index = [index, sampNum - sizeAtom + 1];
    end
    
    GSParam.nSpecialFeat = trainDICParam.isAddLocInfo * 2 + trainDICParam.isAddTimeInfo;
    
    GSParam.index = index;
    GSParam.ncell = length(index);
    [GSParam.R] = bsCreateRMatrix(index, sizeAtom, sampNum);
   
    tmp = zeros(sampNum, sampNum);
    for iCell = 1 : GSParam.ncell
        tmp = tmp + GSParam.R{iCell}' * GSParam.R{iCell};
    end
    GSParam.invTmp = tmp;
    GSParam.invR = inv(gamma * eye(sampNum) + GSParam.invTmp);
%     GSParam.invG = inv(G' * G + gamma * tmp + lambda * eye(sampNum));
    GSParam.omp_G = GSParam.DIC' * GSParam.DIC;
    
    rangeCoef = GSParam.rangeCoef;
    
    switch trainDICParam.normalizationMode
        case {'off', 'none'}
        otherwise
            % normalization
            GSParam.min_values = repmat(rangeCoef(:, 1), 1, GSParam.ncell);
            GSParam.max_values = repmat(rangeCoef(:, 2), 1, GSParam.ncell);
    end
%     GSParam.neiborIndecies = bsGetNeiborIndecies(GSParam.DIC, GSParam.nNeibor);
end
