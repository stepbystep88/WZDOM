function [outputData] = bsPostReBuildByCSR(GInvParam, GSparseInvParam, inputData, options)

    [sampNum, traceNum] = size(inputData);


    GSparseInvParam = bsInitDLSRPkgs(GSparseInvParam, options.gamma, sampNum);

        % tackle the inverse task
    outputData = zeros(sampNum, traceNum);
    dt = GInvParam.dt;
    if GInvParam.isParallel

        pbm = bsInitParforProgress(GInvParam.numWorkers, ...
            traceNum, ...
            'Reconstruct synthetic data progress information', ...
            GInvParam.modelSavePath, ...
            GInvParam.isPrintBySavingFile);

        % parallel computing
        parfor iTrace = 1 : traceNum
            outputData(:, iTrace) = bsHandleOneTrace(GSparseInvParam, inputData(:, iTrace), options, dt);
            bsIncParforProgress(pbm, iTrace, 20000);
        end

        try
            delete(pbm.name);
        catch
        end
    else
        % non-parallel computing 
        for iTrace = 1 : traceNum
%             fprintf('Reconstructing the %d-th trace...\n', iTrace);
            outputData(:, iTrace) = bsHandleOneTrace(GSparseInvParam, inputData(:, iTrace), options, dt);
        end
    end
    
end

function newData = bsHandleOneTrace(GSparseInvParam, realData, options, dt)

    sampNum = size(realData, 1);
    ncell = GSparseInvParam.ncell;
    sizeAtom = GSparseInvParam.sizeAtom;
    patches = zeros(sizeAtom, ncell);
    rangeCoef = options.rangeCoef;
    
    for j = 1 : ncell
        js = GSparseInvParam.index(j);
        patches(:, j) = realData(js : js+sizeAtom-1);
    end
    
    if strcmp(options.mode, 'low_high')
        patches = (patches - rangeCoef(1, 1)) / (rangeCoef(1, 2) - rangeCoef(1, 1));
    end
        
    gammas = omp(GSparseInvParam.D1'*patches, ...
                    GSparseInvParam.omp_G, ...
                    GSparseInvParam.sparsity);
%     gammas = gammas .* GSparseInvParam.C;
    
    new_patches = GSparseInvParam.D2 *  gammas;
    if strcmp(options.mode, 'low_high')
        new_patches = new_patches * (rangeCoef(2, 2) - rangeCoef(2, 1)) + rangeCoef(2, 1); 
    end

    switch GSparseInvParam.reconstructType
        case 'equation'
            avgLog = options.gamma * realData;
            % get reconstructed results by equation
            for j = 1 : ncell

                avgLog = avgLog + GSparseInvParam.R{j}' * new_patches(:, j);
            end

            tmpData = GSparseInvParam.invR * avgLog;
        case 'simpleAvg'
            % get reconstructed results by averaging patches
            avgLog = bsAvgPatches(new_patches, GSparseInvParam.index, sampNum);
            tmpData = avgLog * options.gamma + realData * (1 - options.gamma);
    end
    
    % 合并低频和中低频
    if strcmp(options.mode, 'low_high')
        ft = 1/dt*1000/2;
        newData = bsMixTwoSignal(realData, tmpData, options.lowCut*ft, options.lowCut*ft, dt/1000);
%         bsShowFFTResultsComparison(2, [realData, tmpData, newData], {'反演结果', '高频', '合并'});
    else
        newData = tmpData;
    end
    
end

function GSparseInvParam = bsInitDLSRPkgs(GSparseInvParam, gamma, sampNum)

    validatestring(string(GSparseInvParam.reconstructType), {'equation', 'simpleAvg'});
    
    [sizeAtom, nAtom] = size(GSparseInvParam.DIC);
    sizeAtom = sizeAtom / 2;
    
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
    
    D1 = GSparseInvParam.DIC(1:sizeAtom, :);
    D2 = GSparseInvParam.DIC(sizeAtom+1:end, :);
    
    [D1, D2, C] = bsNormalDIC(D1, D2);
    
    GSparseInvParam.D1 = D1;
    
    GSparseInvParam.omp_G = D1' * D1;
    GSparseInvParam.D2 = D2;
    GSparseInvParam.C = C;
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
  