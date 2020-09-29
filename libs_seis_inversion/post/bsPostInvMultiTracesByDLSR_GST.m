function [data, ys] = bsPostInvMultiTracesByDLSR_GST(GInvParam, neiboors, ds, G, xs, scaleFactors, inIds, crossIds, method)
    % 第二步：大循环
    options = method.options;
    seisOption = GInvParam.seisInvOptions;
    GBOptions = seisOption.GBOptions;
    GBOptions.maxIter = options.innerIter;

	initRegParam = options.initRegParam;
    
    [sampNum, traceNum] = size(xs);
    
    GSParam = method.parampkgs;
    gst_options = GSParam.gst_options;
    
    gamma = method.regParam.gamma;
    lambda = method.regParam.lambda;
    
    xs_org = xs;
    
    if isfield(GSParam, 'nMultipleTrace')
        nBlock = GSParam.nMultipleTrace;
    else
        nBlock = 1;
    end
    
    if isfield(GSParam, 'is3D') && GSParam.is3D
        is3D = true;
        GSParam.nInline = max(inIds) - min(inIds) + 1;
        GSParam.nCrossline = max(crossIds) - min(crossIds) + 1;
    else
        is3D = false;
        GSParam.nInline = [];
        GSParam.nCrossline = [];
    end
    
    
    GSParam = bsInitGSparseParam(GSParam, sampNum, nBlock);
    
    pbm = bsInitParforProgress(GInvParam.numWorkers, ...
            traceNum, ...
            '', ...
            GInvParam.modelSavePath, ...
            GInvParam.isPrintBySavingFile);
        
    if is3D
        
        fprintf('Smooth the seismic data by using GST...\n');
        ds = bsReshapeDataAs3D(ds, GSParam.nInline, GSParam.nCrossline);
%         gst_options.iterNum = 30;
        ds_3D = bsSmoothByGST3D(ds, [], gst_options);
        ds = bsReshapeDataAs2D(ds_3D);
        
        old_d = ds .* repmat(scaleFactors, sampNum-1, 1);
        old_d = [old_d; old_d(end, :)];
        shiftedData = bsPhase90Shift(old_d);     % 相移90°
        
        fprintf('Extract the structure tensor of smoothed seismic data as reference...\n');
        
        % 计算结构张量
        ds_3D = bsReshapeDataAs3D(shiftedData, GSParam.nInline, GSParam.nCrossline);
        S = bsGetStructureTensor3D(ds_3D, gst_options);
        
    else
        fprintf('Smooth the seismic data by using GST...\n');
        
        ds = bsSmoothByGST2D(ds, [], gst_options);
        old_d = ds .* repmat(scaleFactors, sampNum-1, 1);
        old_d = [old_d; old_d(end, :)];
        shiftedData = bsPhase90Shift(old_d);     % 相移90°
        
        fprintf('Extract the structure tensor of smoothed seismic data as reference...\n');
        
        % 计算结构张量
        S = bsGetStructureTensor2D(shiftedData, gst_options);
    end
    
    % 计算G'*G
    sG = G / scaleFactors(1);
    GTG = sG'* sG;
    scaleMat = 1 ./ (repmat(scaleFactors, sampNum, 1) / scaleFactors(1));
    IdentityMat = sparse(eye(sampNum) * (initRegParam + lambda));
    
    %% 利用地震数据反演需要用到的子函数
    function s = bsAx(r)
        s = (GTG * r) .* (scaleMat.^2) + IdentityMat * r;
%         s = s + bsGSTForwardD2D(S, r, gst_alpha);
%         u = nabla(r);
%         Su = Mult(S, u);
% %         Su = u;
% %         uSu = (Su(:,:,1) .* u(:, :, 1) + Su(:,:,2) .* u(:,:,2));
%         s = s + gst_options.tau * div(Su);
    end
    
    function b = bsB(xs_init, xs_org)
        b = (sG' * ds) .* scaleMat + initRegParam * xs_init + lambda * xs_org;
    end
    
    for iter = 1 : options.maxIter
        pbm = bsResetParforProgress(pbm, sprintf("The %d/%d-th iteration: regular inversion.", iter, options.maxIter));

        b = bsB(xs_org, xs);
        xs = bsCGMultiTraces(xs, b, @bsAx, options.innerIter, S, gst_options, is3D, GSParam.nInline, GSParam.nCrossline);
        
        Ips = exp(xs);
    
        pbm = bsResetParforProgress(pbm, sprintf("The %d/%d-th iteration: sparse reconstruction.", iter, options.maxIter));
        
%         newIps = Ips;
        for iTrace = 1 : traceNum
%             avg_xs(:, iTrace) = sparseRebuildOneTrace(GSParam, Ips(:, neiboors{iTrace}));
            tmp = neiboors{iTrace};
            out = bsSparseRebuildOneTrace(GSParam, {Ips(:, tmp)}, gamma, inIds(tmp), crossIds(tmp));
            
            xs(:, iTrace) = log(out);
            bsIncParforProgress(pbm, iTrace, 1000);
        end
    end

    data = exp(xs);
    if ~isempty(gst_options.filter_fcn) && isa(gst_options.filter_fcn, 'function_handle')
        data = gst_options.filter_fcn(data);
    end
    
    ys = ds - G * xs;
    
    bsDeleteParforProgress(pbm);
end
    



