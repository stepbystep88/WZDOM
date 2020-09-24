function [data, ys] = bsPostInvMultiTracesByDLSR_NLM(GInvParam, neiboors, nlm_ps, ds, G, xs, scaleFactors, inIds, crossIds, method)
    % 第二步：大循环
    options = method.options;
    seisOption = GInvParam.seisInvOptions;
    GBOptions = seisOption.GBOptions;
    GBOptions.maxIter = options.innerIter;
        
    mainFunc = seisOption.mainFunc;
	initRegParam = options.initRegParam;
    
    [sampNum, traceNum] = size(xs);
    
    gamma = method.regParam.gamma;
    lambda = method.regParam.lambda;
    xs_org = xs;
    
    GSParam = method.parampkgs;
    nBlock = size(nlm_ps{1}, 1);
    
    GSParam = bsInitGSparseParam(GSParam, sampNum, nBlock);
    
    pbm = bsInitParforProgress(GInvParam.numWorkers, ...
            traceNum, ...
            '', ...
            GInvParam.modelSavePath, ...
            GInvParam.isPrintBySavingFile);
        
    for iter = 1 : options.maxIter
        pbm = bsResetParforProgress(pbm, sprintf("The %d/%d-th iteration: regular inversion.", iter, options.maxIter));
        parfor iTrace = 1 : traceNum
            xs(:, iTrace) = invNormalOneTrace(ds(:, iTrace), G/scaleFactors(iTrace), xs(:, iTrace), ...
                xs_org(:, iTrace), mainFunc, lambda, initRegParam, GBOptions);
            bsIncParforProgress(pbm, iTrace, 200);
        end

        Ips = exp(xs);
    
        pbm = bsResetParforProgress(pbm, sprintf("The %d/%d-th iteration: sparse reconstruction.", iter, options.maxIter));
        for iTrace = 1 : traceNum
%             avg_xs(:, iTrace) = sparseRebuildOneTrace(GSParam, Ips(:, neiboors{iTrace}));
            tmp = neiboors{iTrace};
            out = bsSparseRebuildOneTrace(GSParam, {Ips(:, tmp)}, gamma, inIds(tmp), crossIds(tmp), nlm_ps{iTrace});
            xs(:, iTrace) = log(out);
            bsIncParforProgress(pbm, iTrace, 200);
        end
        
    end

    data = exp(xs);
    ys = ds - G * xs;
end
    

function xOut = invNormalOneTrace(d, G, x, xInit, mainFunc, lambda, initRegParam, GBOptions)
    inputObjFcnPkgs = {
        mainFunc,    struct('A', G, 'B', d),   1; 
        @bsReg1DTKInitModel,    struct('xInit', []), lambda;
        @bsReg1DTKInitModel,    struct('xInit', xInit), initRegParam;
    };
    [xOut, ~, ~, ~] = bsGBSolveByOptions(inputObjFcnPkgs, x, [], [], GBOptions);
end




