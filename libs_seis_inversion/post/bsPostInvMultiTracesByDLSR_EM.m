function [data, ys] = bsPostInvMultiTracesByDLSR_EM(GInvParam, neiboors, ds, G, xs, scaleFactors, inIds, crossIds, method)
    % 第二步：大循环
    options = method.options;
    seisOption = GInvParam.seisInvOptions;
    GBOptions = seisOption.GBOptions;
    GBOptions.maxIter = options.innerIter;
        
    mainFunc = seisOption.mainFunc;
	initRegParam = options.initRegParam;
    
    [sampNum, traceNum] = size(xs);
    
    ys = zeros(sampNum-1, traceNum);
    xs_org = xs;
    
    ds = ds .* repmat(scaleFactors, sampNum-1, 1);
    
    
    
%     ds = bsSmoothByGST2D(ds, [], GST_options);
    
    ds_org = ds;
    
    gamma = method.regParam.gamma;
    lambda = method.regParam.lambda;
    beta = method.regParam.beta;
    
    [GSParam, flag] = bsInitDLSRPkgs(method.parampkgs, sampNum);
    
    if ~isfield(GSParam, 'ratio_to_reconstruction')
        GSParam.ratio_to_reconstruction = 0.1;
    end
    
    if ~isfield(GSParam, 'wellPos')
        GSParam.wellPos = [];
    end
    
    n = traceNum * GSParam.ratio_to_reconstruction;
    seqs = unique([round(linspace(1, traceNum, n)), GSParam.wellPos]);
    seqs = sort(seqs);
    
    [X, Y] = meshgrid(seqs, 1:sampNum-1);
	[Xq,Yq] = meshgrid(1:traceNum, 1:sampNum-1);
    
    [XX, YY] = meshgrid(seqs, 1:sampNum);
	[XXq,YYq] = meshgrid(1:traceNum, 1:sampNum);
    
    pbm = bsInitParforProgress(GInvParam.numWorkers, ...
            traceNum, ...
            '', ...
            GInvParam.modelSavePath, ...
            GInvParam.isPrintBySavingFile);
        
    
    GST_options1 = bsCreateGSTParam(2, 'sigma', 100, 'iterNum', 30);
    GST_options2 = bsCreateGSTParam(2, 'sigma', 10, 'iterNum', 30);
    
    xs = bsSmoothByGST2D(xs, [], GST_options1);
    
    switch flag
        case 2
            GModel = GSParam.model;
            GError = GSParam.error;
            
            for iTrace = seqs
                % 初始化y
                tmp = neiboors{iTrace};
                ys(:, iTrace) = bsSparsePredictOneTrace(GError, {ds_org(:, tmp)}, inIds(tmp), crossIds(tmp));
%                 ds(:, iTrace) = ds_org(:, iTrace) - ys(:, iTrace) * beta;
            end
            
            % 给未预测的道插值
            if GSParam.ratio_to_reconstruction < 1 && traceNum > 1
                ys = interp2(X, Y, ys(:, seqs), Xq, Yq, 'spline');
            end
            ys = bsSmoothByGST2D(ys, ds_org, GST_options2);
            ds = ds_org - ys * beta;
%             ds = ds_org;
            
            for iter = 1 : options.maxIter
                pbm = bsResetParforProgress(pbm, sprintf("The %d-th iteration：regular inversion.", iter));
                parfor iTrace = 1 : traceNum
                    [xs(:, iTrace), ys(:, iTrace)] ...
                        = invNormalOneTrace(...
                        ds(:, iTrace), ...
                        G, ...
                        xs(:, iTrace), ...
                        xs_org(:, iTrace), ...
                        scaleFactors(iTrace), ...
                        lambda, initRegParam, GBOptions);
                    bsIncParforProgress(pbm, iTrace, 1000);
                end

                Ips = exp(xs);
                newIps = Ips;
                
                pbm = bsResetParforProgress(pbm, sprintf("The %d-th iteration：sparse reconstruction.", iter));
                for iTrace = 1 : traceNum
                    tmp = neiboors{iTrace};

                    newIps(:, iTrace) = bsSparseRebuildOneTrace(GModel, {Ips(:, tmp)}, 1, inIds(tmp), crossIds(tmp));
                    bsIncParforProgress(pbm, iTrace, 1000);
                end
                
                xs = log( Ips * (1- gamma) + newIps * gamma );
                
%                 xs = bsSmoothByGST2D(xs, [ds_org; ds_org(end, :)], GST_options1);
                
%                 ys = ds_org - G * log( Ips * (1- beta) + newIps * beta );
                ys = G * xs;
                ys_bak = ys;
                
                for iTrace = seqs
                    % 初始化y
                    tmp = neiboors{iTrace};
                    out = bsSparseRebuildOneTrace(GError, {ys_bak(:, tmp), ds_org(:, tmp)}, 1, inIds(tmp), crossIds(tmp));
                    ys(:, iTrace) = out(:, 1);
                end

                % 给未预测的道插值
                if GSParam.ratio_to_reconstruction < 1 && traceNum > 1
                    ys = interp2(X, Y, ys(:, seqs), Xq, Yq, 'spline');
                end
                
%                 ys = bsFilterProfileData(ys, 0.1, 1);
                
                ys = bsSmoothByGST2D(ys, ds_org, GST_options2);
                ys = ys * beta + ys_bak * (1 - beta);
                ds = ds_org - ys * beta;
            end
        case 1
            nDic = (size(GSParam.DIC, 1) - GSParam.nSpecialFeat) / GSParam.sizeAtom;
            
%             for iTrace = seqs
%                 % 初始化y
%                 tmp = neiboors{iTrace};
%                 inity = bsSparsePredictOneTrace(GSParam, {[ds_org(:, tmp); ds_org(end, tmp)]}, inIds(tmp), crossIds(tmp));
%                 ys(:, iTrace) = inity(1:sampNum-1, 1);
%             end
%             
%             % 给未预测的道插值
%             if GSParam.ratio_to_reconstruction < 1 && traceNum > 1
%                 ys = interp2(X, Y, ys(:, seqs), Xq, Yq, 'spline');
%             end
%             
%             ys = bsSmoothByGST2D(ys, [], GST_options2);
%             ds = ds_org - ys * beta;  
            ds = ds_org;
            
            for iter = 1 : options.maxIter
                pbm = bsResetParforProgress(pbm, sprintf("The %d-th iteration：regular inversion.", iter));
                parfor iTrace = 1 : traceNum
                    [xs(:, iTrace), ys(:, iTrace)] ...
                        = invNormalOneTrace(...
                        ds(:, iTrace), ...
                        G, ...
                        xs(:, iTrace), ...
                        xs_org(:, iTrace), ...
                        scaleFactors(iTrace), ...
                        lambda, initRegParam, GBOptions);
                    bsIncParforProgress(pbm, iTrace, 1000);
                end

                Ips = exp(xs);

                pbm = bsResetParforProgress(pbm, sprintf("The %d-th iteration：sparse reconstruction.", iter));
                ys_bak = ys;
                for iTrace = 1:traceNum
                    tmp = neiboors{iTrace};
                    
                    switch nDic
                        case 3
                            input = {Ips(:, tmp), [ys_bak(:, tmp); ys_bak(end, tmp)], [ds_org(:, tmp); ds_org(end, tmp)]};
                        case 2
                            input = {Ips(:, tmp), [ys_bak(:, tmp); ys_bak(end, tmp)]};
                    end
                
                    out = bsSparseRebuildOneTrace(GSParam, input, 1, inIds(tmp), crossIds(tmp));
                    
                    xs(:, iTrace) = out(:, 1); 
                    ys(:, iTrace) = out(1:sampNum-1, 2);
                    
                    bsIncParforProgress(pbm, iTrace, 1000);
                end
                
            
                % 给未预测的道插值
                if GSParam.ratio_to_reconstruction < 1 && traceNum > 1
                    ys = interp2(X, Y, ys(:, seqs), Xq, Yq, 'spline');
%                     xs = interp2(XX, YY, xs(:, seqs), XXq, YYq, 'spline');
                end
                
%                 xs = bsSmoothByGST2D(xs, [], GST_options1);
                xs = (xs  * gamma + (1 - gamma) * Ips );
                xs = log(xs);
                
                ys = ys * beta + ys_bak * (1 - beta);
                ys = bsSmoothByGST2D(ys, [], GST_options2);
                ds = ds_org - ys * beta;
            end
    end
    

    data = exp(xs);
    ys = ds_org - ds;
end
    

function [xOut, y] = invNormalOneTrace(d, G, x, xInit, scale, lambda, initRegParam, GBOptions)
    inputObjFcnPkgs = {
        @bsLinearTwoNorm,    struct('A', G/scale, 'B', d/scale),   1; 
        @bsReg1DTKInitModel,    struct('xInit', []), lambda;
        @bsReg1DTKInitModel,    struct('xInit', xInit), initRegParam;
    };

    [xOut, ~, ~, ~] = bsGBSolveByOptions(inputObjFcnPkgs, x, [], [], GBOptions);
    y = d - (G * xOut);
%     y = [y; y(end)];
end

function [GSParam, flag] = bsInitDLSRPkgs(GSParam, sampNum)
    
    if isfield(GSParam, 'error') &&  isfield(GSParam, 'model')
        flag = 2;
        
        if isfield(GSParam.model, 'omp_G')
            return;
        end
        
        if isfield(GSParam, 'nMultipleTrace')
            nBlock = GSParam.nMultipleTrace;
        else
            nBlock = 1;
        end
    
        [GSParam.model] = bsInitGSparseParam(GSParam.model, sampNum, nBlock);
        [GSParam.error] = bsInitGSparseParam(GSParam.error, sampNum - 1, nBlock, [], 1);
        
    else
        flag = 1;
        
        if isfield(GSParam, 'omp_G')
            return;
        end
        
        if isfield(GSParam, 'nMultipleTrace')
            nBlock = GSParam.nMultipleTrace;
        else
            nBlock = 1;
        end
    
        nDic = floor(size(GSParam.DIC, 1) / GSParam.trainDICParam.sizeAtom);
        switch nDic
            case 2
                [GSParam] = bsInitGSparseParam(GSParam, sampNum, nBlock, [], 1);
            case 3
                [GSParam] = bsInitGSparseParam(GSParam, sampNum, nBlock, 3, 2);
        end
        
    end
end
