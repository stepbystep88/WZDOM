function [vp, vs, rho] = bsPreInvMultiTracesByCSR(GInvParam, neiboors, ds, G, xs, scaleFactors, lsdCoef, inIds, crossIds, method)
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
    
%     if isfield(GSParam, 'nMultipleTrace')
%         nBlock = GSParam.nMultipleTrace;
%     else
%         nBlock = 1;
%     end
    nBlock = length(neiboors{1});
    
    GSParam = bsInitGSparseParam(GSParam, sampNum/3, nBlock);
    nDic = (size(GSParam.DIC, 1) - GSParam.nSpecialFeat) / GSParam.sizeAtom;
    mode = GInvParam.mode;
%     lsdCoef = GInvParam.lsdCoef;
    
    pbm = bsInitParforProgress(GInvParam.numWorkers, ...
            traceNum, ...
            '', ...
            GInvParam.modelSavePath, ...
            GInvParam.isPrintBySavingFile);
        
    for iter = 1 : options.maxIter
        pbm.title = sprintf("The %d/%d-th iteration: regular inversion] ", iter, options.maxIter);
        parfor iTrace = 1 : traceNum
            xs(:, iTrace) = invNormalOneTrace(ds(:, iTrace), G/scaleFactors(iTrace), xs(:, iTrace), ...
                xs_org(:, iTrace), mainFunc, lambda, initRegParam, GBOptions);
            
            bsIncParforProgress(pbm, iTrace, 200);
            
        end

        [vp, vs, rho] = bsPreRecoverElasticParam(xs, mode, lsdCoef);
        vp_vs = vp ./ vs;
        
        pbm.title = sprintf("The %d/%d-th iteration: sparse reconstruction ", iter, options.maxIter);
%         fprintf("第%d次迭代：稀疏重构\n", iter);
        for iTrace = 1 : traceNum
%             avg_xs(:, iTrace) = sparseRebuildOneTrace(GSParam, Ips(:, neiboors{iTrace}));
            tmp = neiboors{iTrace};
            
            switch nDic
                case 3
                    [out, ~] = bsSparseRebuildOneTrace(GSParam, {vp(:, tmp), vs(:, tmp), rho(:, tmp)}, gamma, inIds(tmp), crossIds(tmp));
                    newData = [out(:, 1), out(:, 1 : 3)];
                case 4
                    [out, ~] = bsSparseRebuildOneTrace(GSParam, {vp(:, tmp), vs(:, tmp), rho(:, tmp), vp_vs(:, tmp)}, gamma, inIds(tmp), crossIds(tmp));
                    newData = [out(:, 1), out(:, 1), out(:, 1)./out(:, 4), out(:, 3)];
            end
            
            
        
            xs(:, iTrace) = bsPreBuildModelParam(newData, mode, lsdCoef{iTrace});
            
            bsIncParforProgress(pbm, iTrace, 200);
        end

    end

    [vp, vs, rho] = bsPreRecoverElasticParam(xs, mode, lsdCoef);
%     ys = ds - G * xs;
end
    

function xOut = invNormalOneTrace(d, G, x, xInit, mainFunc, lambda, initRegParam, GBOptions)
    inputObjFcnPkgs = {
        mainFunc,    struct('A', G, 'B', d),   1; 
        @bsReg1DTKInitModel,    struct('xInit', []), lambda;
        @bsReg1DTKInitModel,    struct('xInit', xInit), initRegParam;
    };
    [xOut, ~, ~, ~] = bsGBSolveByOptions(inputObjFcnPkgs, x, [], [], GBOptions);
end




