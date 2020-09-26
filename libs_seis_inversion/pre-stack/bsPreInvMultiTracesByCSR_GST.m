function [vp, vs, rho] = bsPreInvMultiTracesByCSR_GST(GInvParam, neiboors, ds, G, xs, scaleFactors, lsdCoef, shiftedData, inIds, crossIds, method)
    % 第二步：大循环
    options = method.options;
    seisOption = GInvParam.seisInvOptions;
    GBOptions = seisOption.GBOptions;
    GBOptions.maxIter = options.innerIter;
     
    mainFunc = seisOption.mainFunc;
	initRegParam = options.initRegParam;
    
    [sampNum, traceNum] = size(xs);
    xs_org = xs;
    
    GSParam = method.parampkgs;
    
    if isfield(GSParam, 'nMultipleTrace')
        nBlock = GSParam.nMultipleTrace;
    else
        nBlock = 1;
    end
    
    GSParam = bsInitGSparseParam(GSParam, sampNum/3, nBlock);
    nDic = (size(GSParam.DIC, 1) - GSParam.nSpecialFeat) / GSParam.sizeAtom;
    mode = GInvParam.mode;
%     lsdCoef = GInvParam.lsdCoef;
    
    gamma = method.regParam.gamma;
    lambda = method.regParam.lambda;
    
    if isfield(GSParam, 'is3D') && GSParam.is3D
        is3D = true;
        GSParam.nInline = max(inIds) - min(inIds) + 1;
        GSParam.nCrossline = max(crossIds) - min(crossIds) + 1;
    else
        is3D = false;
        GSParam.nInline = [];
        GSParam.nCrossline = [];
    end
    
    pbm = bsInitParforProgress(GInvParam.numWorkers, ...
            traceNum, ...
            '', ...
            GInvParam.modelSavePath, ...
            GInvParam.isPrintBySavingFile);
        
    %% 结构张量相关
    if isfield(GSParam, 'gst_options')
        gst_options = GSParam.gst_options;
    else
        if is3D
            gst_options = bsCreateGSTParam(3);
        else
            gst_options = bsCreateGSTParam(2);
        end
        
    end
    
    
    % 计算G' * G
%     GTGs = cell(1, traceNum);
%     parfor iTrace = 1 : traceNum
%         GTGs{iTrace} = Gs{iTrace}' * Gs{iTrace};
%     end
%     
%     IdentityMat = sparse(eye(sampNum) * (initRegParam + lambda));
    
    % 计算结构张量信息
    if is3D
        gst_options.iterNum = 50;
        shiftedData = bsSmoothByGST3D(bsReshapeDataAs3D(shiftedData, GSParam.nInline, GSParam.nCrossline), [], gst_options);
        
        [S, ~, blur] = bsGetStructureTensor3D(shiftedData, gst_options);
    else
        shiftedData = bsSmoothByGST2D(shiftedData, [], gst_options);
        [S, ~, ~, ~, ~,~, blur] = bsGetStructureTensor2D(shiftedData, gst_options);
    end
    
    
     %% 利用地震数据反演需要用到的子函数
%     function s = bsAx(r)
%         s = IdentityMat * r;
% %         parfor i = 1 : traceNum
% %             s(:, i) = s(:, i) + GTGs{i} * r(:, i);
% %         end
%         s = s +  GTGs{1} * r .* (repmat(scaleFactors, sampNum, 1) / scaleFactors(1));
%         s(1:sampNum/3, :) = s(1:sampNum/3, :) + bsGSTForwardD2D(D, r(1:sampNum/3, :), gst_alpha);
%     end
%     
%     function b = bsB(xs_init, xs_org)
%         b = initRegParam * xs_init + lambda * xs_org;
%         for i = 1 : traceNum
%             b(:, i) = b(:, i) + Gs{i}' * ds(:, i);
%         end
%     end
    gst_options.iterNum = options.innerIter;
    for iter = 1 : options.maxIter
        pbm.title = sprintf("The %d/%d-th iteration: regular inversion.", iter, options.maxIter);
%         b = bsB(xs_org, xs);
%         xs = bsCGMultiTraces(xs, b, @bsAx, options.innerIter);
%         G=Gs{1}*scaleFactors(1);
        
        parfor iTrace = 1 : traceNum
            xs(:, iTrace) = invNormalOneTrace(ds(:, iTrace), G/scaleFactors(iTrace), xs(:, iTrace), ...
                xs_org(:, iTrace), mainFunc, lambda, initRegParam, GBOptions);
            
            bsIncParforProgress(pbm, iTrace, 200);
        end
        
        [vp, vs, rho] = bsPreRecoverElasticParam(xs, mode, lsdCoef);
        
        % 利用gst光滑
        fprintf("[The %d/%d-th iteration: smooth results.]\n", iter, options.maxIter);
        
%         if iter == options.maxIter
        if is3D
%             vp = bsSmoothByGST3D(vp, [], gst_options, S, blur);
%             vs = bsSmoothByGST3D(vs, [], gst_options, S, blur);
%             rho = bsSmoothByGST3D(rho, [], gst_options, S, blur);
            vp = bsReshapeDataAs2D(bsSmoothByGST3D(bsReshapeDataAs3D(vp, GSParam.nInline, GSParam.nCrossline), [], gst_options, S, blur));
            vs = bsReshapeDataAs2D(bsSmoothByGST3D(bsReshapeDataAs3D(vs, GSParam.nInline, GSParam.nCrossline), [], gst_options, S, blur));
            rho = bsReshapeDataAs2D(bsSmoothByGST3D(bsReshapeDataAs3D(rho, GSParam.nInline, GSParam.nCrossline), [], gst_options, S, blur));
        else
            vp = bsSmoothByGST2D(vp, [], gst_options, S, blur);
            vs = bsSmoothByGST2D(vs, [], gst_options, S, blur);
            rho = bsSmoothByGST2D(rho, [], gst_options, S, blur);
        end
%         end
        
        vp_vs = vp ./ vs;
        
        pbm.title = sprintf("The %d/%d-th iteration: sparse reconstruction ", iter, options.maxIter);
%         fprintf("第%d次迭代：稀疏重构\n", iter);
        if nBlock == 1
            parfor iTrace = 1 : traceNum
                switch nDic
                    case 3
                        [out, ~] = bsSparseRebuildOneTrace(GSParam, {vp(:, iTrace), vs(:, iTrace), rho(:, iTrace)}, gamma, inIds(iTrace), crossIds(iTrace));
                    case 4
                        [out, ~] = bsSparseRebuildOneTrace(GSParam, {vp(:, iTrace), vs(:, iTrace), rho(:, iTrace), vp_vs(:, iTrace)}, gamma, inIds(iTrace), crossIds(iTrace));
                end

                xs(:, iTrace) = bsPreBuildModelParam([out(:, 1), out(:, 1 : 3)], mode, lsdCoef{iTrace});

                bsIncParforProgress(pbm, iTrace, 200);
            end
            
        else
            for iTrace = 1 : traceNum
    %             avg_xs(:, iTrace) = sparseRebuildOneTrace(GSParam, Ips(:, neiboors{iTrace}));
                tmp = neiboors{iTrace};

                switch nDic
                    case 3
                        [out2, ~] = bsSparseRebuildOneTrace(GSParam, {vp(:, tmp), vs(:, tmp), rho(:, tmp)}, gamma, inIds(tmp), crossIds(tmp));
                    case 4
                        [out2, ~] = bsSparseRebuildOneTrace(GSParam, {vp(:, tmp), vs(:, tmp), rho(:, tmp), vp_vs(:, tmp)}, gamma, inIds(tmp), crossIds(tmp));
                end

                xs(:, iTrace) = bsPreBuildModelParam([out2(:, 1), out2(:, 1 : 3)], mode, lsdCoef{iTrace});

                bsIncParforProgress(pbm, iTrace, 200);
            end
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




