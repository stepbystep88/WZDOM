function [data, ys] = bsPostInvMultiTracesByDLSR(GInvParam, neiboors, ds, G, xs, scaleFactors, inIds, crossIds, method)
    % 第二步：大循环
    options = method.options;
    seisOption = GInvParam.seisInvOptions;
    GBOptions = seisOption.GBOptions;
    GBOptions.maxIter = options.innerIter;
        
    mainFunc = seisOption.mainFunc;
	initRegParam = options.initRegParam;
    
    [sampNum, traceNum] = size(xs);
    
    avg_xs = zeros(sampNum, traceNum);
    gamma = method.regParam.gamma;
    lambda = method.regParam.lambda;
    xs_org = xs;
    
    GSParam = method.parampkgs;
    
    if isfield(GSParam, 'nMultipleTrace')
        nBlock = GSParam.nMultipleTrace;
    else
        nBlock = 1;
    end
    
    GSParam = bsInitGSparseParam(GSParam, sampNum, nBlock);
    
    for iter = 1 : options.maxIter
        fprintf("第%d次迭代：常规反演\n", iter);
        parfor iTrace = 1 : traceNum
            xs(:, iTrace) = invNormalOneTrace(ds(:, iTrace), G/scaleFactors(iTrace), xs(:, iTrace), ...
                xs_org(:, iTrace), mainFunc, lambda, initRegParam, GBOptions);
        end

        Ips = exp(xs);

        fprintf("第%d次迭代：稀疏重构\n", iter);
        for iTrace = 1 : traceNum
%             avg_xs(:, iTrace) = sparseRebuildOneTrace(GSParam, Ips(:, neiboors{iTrace}));
            tmp = neiboors{iTrace};
            out = bsSparseRebuildOneTrace(GSParam, {Ips(:, tmp)}, gamma, inIds(tmp), crossIds(tmp));
            xs(:, iTrace) = log(out);
        end

%         fprintf("第%d次迭代：反演结果合并\n", iter);
%         parfor iTrace = 1 : traceNum
%             xNew = Ips(:, iTrace) * (1 - gamma) + avg_xs(:, iTrace) * gamma;
%             xs(:, iTrace) = log(xNew);
%         end

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




