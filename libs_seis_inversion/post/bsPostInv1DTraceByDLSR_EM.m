function [x, fval, exitFlag, output] = bsPostInv1DTraceByDLSR_EM(d, G, xInit, Lb, Ub, regParam, parampkgs, options)
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
%     mainData.y = zeros(size(d));
    
    sampNum = length(xInit);
    
    inline = options.inline;
    crossline = options.crossline;
    scaleFactor = options.scaleFactor;
    
    % re-organize the input objective function pakages
    inputObjFcnPkgs = {
        @bsLinearTwoNorm,       mainData,   1; 
        @bsReg1DTKInitModel,    struct('xInit', xInit), options.initRegParam;
        @bsReg1DTKInitModel,    struct('xInit', xInit), options.initRegParam;
    };
    
    % 判定含有lambda、beta和gamma参数
    assert(isfield(regParam, 'lambda') && isfield(regParam, 'beta') && isfield(regParam, 'gamma'))
    
           
    GBOptions = options.GBOptions;
    inputObjFcnPkgs{2, 3} = regParam.lambda;
    GBOptions.maxIter = options.innerIter;
    
    % create packages for sparse inversion 
    [GSParam, flag] = bsInitDLSRPkgs(parampkgs, sampNum);
    midX = [];
    midF = [];
    
    gamma = regParam.gamma;
    beta = regParam.beta;
    
    
    
    %% two methods
    switch flag
        case 1
            % 初始化y
%             inity = bsSparsePredictOneTrace(GSParam, {[d; d(end)] * scaleFactor}, inline, crossline);
%             inputObjFcnPkgs{1, 2}.B = d - inity(1:sampNum-1, 1) / scaleFactor * beta;
            
            nDic = (size(GSParam.DIC, 1) - GSParam.nSpecialFeat) / GSParam.sizeAtom;
            
            for i = 1 : options.maxIter
        
                % change the current initial guess
                inputObjFcnPkgs{2, 2} = [];
                [xOut, fval, exitFlag, output_] = bsGBSolveByOptions(inputObjFcnPkgs, xInit, Lb, Ub, GBOptions);

                if GBOptions.isSaveMiddleRes
                    midX = [midX, output_.midResults.x];
                    midF = [midF, output_.midResults.f];
                end

                es = d - G*xOut;
                switch nDic
                    case 3
                        input = {exp(xOut), [es; es(end)] * scaleFactor, [d; d(end)] * scaleFactor};
                    case 2
                        input = {exp(xOut), [es; es(end)] * scaleFactor};
                end
                out = bsSparseRebuildOneTrace(GSParam, input, 1, inline, crossline);

                xNew = log( out(:, 1) * gamma + (1 - gamma) * input{1} );
                inputObjFcnPkgs{1, 2}.B = d - (out(1:sampNum-1, 2) * beta + es * (1 - beta)) / scaleFactor;

                % projection
                xInit = bsPFunc(xNew, output_.options.bounds);
            end
            isSparseRebuild = GSParam.isSparseRebuild;
            
        case 2
            inity = bsSparsePredictOneTrace(GSParam.error, {d * scaleFactor}, inline, crossline);
            inputObjFcnPkgs{1, 2}.B = d - inity / scaleFactor * beta;
            
            for i = 1 : options.maxIter
        
                % change the current initial guess
                inputObjFcnPkgs{2, 2} = [];
                [xOut, fval, exitFlag, output_] = bsGBSolveByOptions(inputObjFcnPkgs, xInit, Lb, Ub, GBOptions);

                if GBOptions.isSaveMiddleRes
                    midX = [midX, output_.midResults.x];
                    midF = [midF, output_.midResults.f];
                end

                % 接下来的代码时对模型m进行稀疏重构
                xNew = bsSparseRebuildOneTrace(GSParam.model, {exp(xOut)}, gamma, inline, crossline);
                xNew = log(xNew);
                
                % 接下来的代码时对噪声y进行稀疏重构
                yNew = bsSparseRebuildOneTrace(GSParam.error, {(d - G * xNew) * scaleFactor, d * scaleFactor}, beta, inline, crossline);
                inputObjFcnPkgs{1, 2}.B = d - yNew(:, 1) / scaleFactor;
                
                % projection
                xInit = bsPFunc(xNew, output_.options.bounds);
            end
            isSparseRebuild = GSParam.model.isSparseRebuild;
    end
    
    %%
    switch isSparseRebuild
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



% function [f, g] = bsLinearTwoNormWithError(x, data)
%     
%     z = data.A * x - (data.B - data.y);
% 
%     f = sum(z.^2, 1);
%     g = 2*(data.A' * z) ;
% end


function [GSParam, flag] = bsInitDLSRPkgs(GSParam, sampNum)
    
    if isfield(GSParam, 'error') &&  isfield(GSParam, 'model')
        flag = 2;
        if isfield(GSParam.model, 'omp_G')
            return;
        end
        [GSParam.model] = bsInitGSparseParam(GSParam.model, sampNum, 1);
        [GSParam.error] = bsInitGSparseParam(GSParam.error, sampNum - 1, 1, [], 1);
 
    else
        flag = 1;
        if isfield(GSParam, 'omp_G')
            return;
        end
        
        nDic = floor(size(GSParam.DIC, 1) / GSParam.trainDICParam.sizeAtom);
        switch nDic
            case 2
                [GSParam] = bsInitGSparseParam(GSParam, sampNum, 1, [], 1);
            case 3
                [GSParam] = bsInitGSparseParam(GSParam, sampNum, 1, 3, 2);
        end
        
    end
end
