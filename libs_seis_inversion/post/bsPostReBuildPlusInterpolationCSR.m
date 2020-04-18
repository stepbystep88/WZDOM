function [outputData, highData, gamma_vals, gamma_locs] = ...
    bsPostReBuildPlusInterpolationCSR(GInvParam, GSParam, inputData, inIds, crossIds, options)

    [sampNum, traceNum] = size(inputData);


    GSParam = bsInitDLSRPkgs(GSParam, options.gamma, sampNum);

        % tackle the inverse task
    outputData = zeros(sampNum, traceNum);
    highData = zeros(sampNum, traceNum);
    gamma_vals = zeros(GSParam.sparsity * GSParam.ncell, traceNum);
    gamma_locs = zeros(GSParam.sparsity * GSParam.ncell, traceNum);
    
    dt = GInvParam.dt;

    pbm = bsInitParforProgress(GInvParam.numWorkers, ...
        traceNum, ...
        'Rebuid data progress information', ...
        GInvParam.modelSavePath, ...
        GInvParam.isPrintBySavingFile);

    n = traceNum * options.ratio_to_reconstruction;
    seqs = round(linspace(1, traceNum, n));
    
    for iTrace = seqs
        [highData(:, iTrace), gamma_vals(:, iTrace), gamma_locs(:, iTrace)] ...
            = bsCalcHighFreqOfOneTrace(GSParam, inputData(:, iTrace), inIds(iTrace), crossIds(iTrace), options);
        bsIncParforProgress(pbm, iTrace, 10000);
    end

    % 给未高分辨率的道插值
    if options.ratio_to_reconstruction < 1 && traceNum > 1
        [X, Y] = meshgrid(seqs, 1:sampNum);
        [Xq,Yq] = meshgrid(1:traceNum, 1:sampNum);

        highData = interp2(X, Y, highData(:, seqs), Xq, Yq, 'spline');
    end
    
    % parallel computing
    for iTrace = 1 : traceNum
        outputData(:, iTrace) = bsHandleOneTrace(inputData(:, iTrace), highData(:, iTrace), options, dt);
        bsIncParforProgress(pbm, iTrace, 10000);
    end

    bsDeleteParforProgress(pbm);
        
end

function [avgLog, gamma_vals, gamma_locs] = bsCalcHighFreqOfOneTrace(GSParam, realData, x, y, options)
    ncell = GSParam.ncell;
    sizeAtom = GSParam.sizeAtom;
    patches = zeros(sizeAtom, ncell);
    rangeCoef = GSParam.rangeCoef;
    
    trainDICParam = GSParam.trainDICParam;
    
    for j = 1 : ncell
        js = GSParam.index(j);
        patches(:, j) = realData(js : js+sizeAtom-1);
    end
    
    if trainDICParam.isAddLocInfo && trainDICParam.isAddTimeInfo
        patches = [ones(1, ncell) * x; ones(1, ncell) * y; 1 : ncell; patches];
    elseif trainDICParam.isAddLocInfo
        patches = [ones(1, ncell) * x; ones(1, ncell) * y; patches];
    elseif trainDICParam.isAddTimeInfo
        patches = [1 : ncell; patches];
    end
    
    switch options.mode
    case {'low_high', 'seismic_high', 'residual'}
    
        switch trainDICParam.normalizationMode
            case 'feat_max_min'
                patches = (patches - GSParam.low_min_values) ./ (GSParam.low_max_values - GSParam.low_min_values);
            case 'whole_data_max_min'
                patches = (patches - rangeCoef(1, 1)) / (rangeCoef(1, 2) - rangeCoef(1, 1));
%             case 'patch_mean_norm'
%                 mean_vlaues = norm(patches, 1);
%                 mean_vlaues = repmat(mean_vlaues, sizeAtom, 1);
%                 patches = (patches - min_values) ./ (max_values - min_values);
        end
        
    end
        
    if strcmp(trainDICParam.feature_reduction, 'high_resolution')
        patches = GSParam.output.B' * patches;
    end
    
    gammas = omp(GSParam.D1'*patches, ...
                    GSParam.omp_G, ...
                    GSParam.sparsity);
            
    [gamma_vals, gamma_locs] = bsGetNonZeroElements(gammas, GSParam.sparsity);
    
    new_patches = GSParam.D2 *  gammas;
    switch options.mode
    case {'low_high', 'seismic_high', 'residual'}
        
        switch trainDICParam.normalizationMode
            case 'feat_max_min'
                new_patches = new_patches .* (GSParam.high_max_values - GSParam.high_min_values) + GSParam.high_min_values; 
            case 'whole_data_max_min'
                new_patches = new_patches .* (rangeCoef(2, 2) - rangeCoef(2, 1)) + rangeCoef(2, 1); 
        end
        
    end

    avgLog = bsAvgPatches(new_patches, GSParam.index, size(realData, 1));

end

function newData = bsHandleOneTrace(realData, avgData, options, dt)

    
    % 合并低频和中低频
    switch options.mode
        case {'low_high', 'seismic_high'}

            ft = 1/dt*1000/2;
            newData = bsMixTwoSignal(realData, avgData, options.lowCut*ft, options.lowCut*ft, dt/1000);
%         bsShowFFTResultsComparison(1, [realData, avgData, newData], {'反演结果', '高频', '合并'});
        case {'full_freq'}
            newData = avgData;
        case 'residual'
            newData = avgData + realData;
    end
    
end



function GSParam = bsInitDLSRPkgs(GSParam, gamma, sampNum)

    validatestring(string(GSParam.reconstructType), {'equation', 'simpleAvg'});
    
    sizeAtom = GSParam.trainDICParam.sizeAtom;
    nAtom = GSParam.trainDICParam.nAtom;

    
    GSParam.sizeAtom = sizeAtom;
    GSParam.nAtom = nAtom;
    GSParam.nrepeat = sizeAtom - GSParam.stride;
    
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
    GSParam.invR = inv(gamma * eye(sampNum) + GSParam.invTmp);
    
    % 低分辨率patch可能有时间和地址信息
    n1 = size(GSParam.DIC, 1) - sizeAtom;
%     nSpecialFeat = trainDICParam.isAddLocInfo * 2 + trainDICParam.isAddTimeInfo;
    
    D1 = GSParam.DIC(1:n1, :);
    D2 = GSParam.DIC(n1+1:end, :);
    
    [D1, D2, C] = bsNormalDIC(D1, D2);
    
    GSParam.D1 = D1;
    
    GSParam.omp_G = D1' * D1;
    GSParam.D2 = D2;
    GSParam.C = C;
    
%     GSParam.neiborIndecies = bsGetNeiborIndecies(D1, GSParam.nNeibor);
    rangeCoef = GSParam.rangeCoef;
    n1 = size(rangeCoef, 1) - sizeAtom;
    
    GSParam.low_min_values = repmat(rangeCoef(1:n1, 1), 1, GSParam.ncell);
    GSParam.low_max_values = repmat(rangeCoef(1:n1, 2), 1, GSParam.ncell);
                
    GSParam.high_min_values = repmat(rangeCoef(n1+1:end, 1), 1, GSParam.ncell);
    GSParam.high_max_values = repmat(rangeCoef(n1+1:end, 2), 1, GSParam.ncell);
end

function [D1, D2, C] = bsNormalDIC(D1, D2)
    C = zeros(size(D1, 2), 1);
    
    for i = 1 : size(D1, 2)
        normCoef = norm(D1(:, i));
        D1(:, i) = D1(:, i) / normCoef;
        D2(:, i) = D2(:, i) / normCoef;
        C(i) = normCoef;
    end
end
  