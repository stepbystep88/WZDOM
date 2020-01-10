function [inIds, crossIds, GInvParam, dstFileNames, segyInfo, type] = bsPostBuildSynData(GInvParam, GSparseInvParam, timeLine, varargin)
%% inverse multiple traces
% Programmed by: Bin She (Email: bin.stepbystep@gmail.com)
% Programming dates: Nov 2019
% -------------------------------------------------------------------------
%
% Input 
% GInvParam     all information of the inverse task
% inIds      	inline ids to be inverted
% crossIds      crossline ids to be inverted
% timeLine      horizon information
% methods       the methods to solve the inverse task, each method is a
% struct depicting the details of the methods. For example, method.load
% indicates whether load the result from mat or segy directly.
% 
% Output
% invVals       inverted results, a cell includ
% -------------------------------------------------------------------------
     p = inputParser;
    
    [rangeInline, rangeCrossline, fileName, segyInfo] ...
        = bsGetWorkAreaRangeByParam(GInvParam);
    
    addParameter(p, 'title', '');
    addParameter(p, 'dstPath', sprintf('%s/seismic/', GInvParam.modelSavePath));
    addParameter(p, 'rangeInline', rangeInline);
    addParameter(p, 'rangeCrossline', rangeCrossline);
    addParameter(p, 'inIds', []);
    addParameter(p, 'crossIds', []);
    addParameter(p, 'expandNum', 0);
    addParameter(p, 'isRebuild', 0);
    addParameter(p, 'gamma', 0.2);
    
    p.parse(varargin{:});  
    options = p.Results;
    
    if isempty(options.inIds) || isempty(options.crossIds)
        [inIds, crossIds] = bsGetCDPsByRange(options.rangeInline, options.rangeCrossline);
    else
        inIds = options.inIds;
        crossIds = options.crossIds;
        
        options.rangeInline = [min(inIds), max(inIds)];
        options.rangeCrossline = [min(crossIds), max(crossIds)];
    end
    
    assert(length(inIds) == length(crossIds), 'The length of inline ids and crossline ids must be the same.');
    
    % create folder to save the intermediate results
    mkdir(options.dstPath);
    type = {'synthetic', 'error_model'};
    dstFileNames{1} = bsGetDstFileName(GSparseInvParam, options, type{1});
    if exist(dstFileNames{1}, 'file') && options.isRebuild == 0
        return;
    else
        GInvParam.upNum = GInvParam.upNum + options.expandNum;
        GInvParam.downNum = GInvParam.downNum + options.expandNum;
        
        traceNum = length(inIds);
        sampNum = GInvParam.upNum + GInvParam.downNum;
    
        % horion of the whole volume
        usedTimeLine = timeLine{GInvParam.usedTimeLineId};
    
        % horizon of given traces
        horizonTimes = bsGetHorizonTime(usedTimeLine, inIds, crossIds, GInvParam.isParallel, GInvParam.numWorkers);

        startTimes = horizonTimes - GInvParam.dt * GInvParam.upNum;

        postSeisData = bsReadTracesByIds(...
            GInvParam.postSeisData.fileName, ...
            GInvParam.postSeisData.segyInfo, ...
            inIds, ...
            crossIds, ...
            startTimes, ...
            sampNum,...
            GInvParam.dt);
    
        GSparseInvParam = bsInitDLSRPkgs(GSparseInvParam, options.gamma, sampNum);
        
        errorData = bsHandleAllTraces();
        synData = postSeisData - errorData;
        
        res.inIds = inIds;
        res.crossIds = crossIds;
        res.horizon = horizonTimes;
        res.upNum = GInvParam.upNum;
        res.dt = GInvParam.dt;
        bsWriteInvResultIntoSegyFile(res, synData, fileName, segyInfo, dstFileNames{1});
        
        dstFileNames{2} = bsGetDstFileName(GSparseInvParam, options, type{2});
        
        bsWriteInvResultIntoSegyFile(res, errorData, fileName, segyInfo, dstFileNames{2});
        
        
    end
    
    function data = bsHandleAllTraces()
        % tackle the inverse task
        data = zeros(sampNum, traceNum);
        gamma = options.gamma;
        
        if GInvParam.isParallel
            
            pbm = bsInitParforProgress(GInvParam.numWorkers, ...
                traceNum, ...
                'Reconstruct synthetic data progress information', ...
                GInvParam.modelSavePath, ...
                GInvParam.isPrintBySavingFile);
            
            % parallel computing
            parfor iTrace = 1 : traceNum
                data(:, iTrace) = bsHandleOneTrace(GInvParam, GSparseInvParam, postSeisData(:, iTrace), gamma);
                bsIncParforProgress(pbm, iTrace, 10000);
            end
            

        else
            % non-parallel computing 
            for iTrace = 1 : traceNum
                fprintf('Reconstructing the trace of inline=%d and crossline=%d...\n', inIds(iTrace), crossIds(iTrace));
                data(:, iTrace) = bsHandleOneTrace(GInvParam, GSparseInvParam, postSeisData(:, iTrace), gamma);
            end
        end
    end
    
end

function dstFileName = bsGetDstFileName(GSparseInvParam, options, type)
    dstFileName = sprintf('%s/%s_%s_gamma=%.2f_sparsity_%d_inline_[%d_%d]_crossline_[%d_%d].sgy', ...
        options.dstPath, type, options.title, ...
        options.gamma, GSparseInvParam.sparsity, ...
        options.rangeInline(1), options.rangeInline(end), ...
        options.rangeCrossline(1), options.rangeCrossline(2));
end


function errorData = bsHandleOneTrace(GInvParam, GSparseInvParam, realData, gamma)

    sampNum = GInvParam.upNum + GInvParam.downNum;
    ncell = GSparseInvParam.ncell;
    sizeAtom = GSparseInvParam.sizeAtom;
    patches = zeros(sizeAtom, ncell);
    
    for j = 1 : ncell
        js = GSparseInvParam.index(j);
        patches(:, j) = realData(js : js+sizeAtom-1);
    end
    
    gammas = omp(GSparseInvParam.D1'*patches, ...
                    GSparseInvParam.omp_G, ...
                    GSparseInvParam.sparsity);
    new_patches = GSparseInvParam.D2 *  gammas;
    

    switch GSparseInvParam.reconstructType
        case 'equation'
            avgLog = gamma * realData;
            % get reconstructed results by equation
            for j = 1 : ncell

                avgLog = avgLog + GSparseInvParam.R{j}' * new_patches(:, j);
            end

            errorData = GSparseInvParam.invR * avgLog;
        case 'simpleAvg'
            % get reconstructed results by averaging patches
            avgLog = bsAvgPatches(new_patches, GSparseInvParam.index, sampNum);
%             errorData = avgLog * gamma + realData * (1 - gamma);
            errorData = avgLog * gamma;
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
    
    [D1, D2] = bsNormalDIC(D1, D2);
    
    GSparseInvParam.D1 = D1;
    
    GSparseInvParam.omp_G = D1' * D1;
    GSparseInvParam.D2 = D2;
end

function [D1, D2] = bsNormalDIC(D1, D2)
    for i = 1 : size(D1, 2)
        normCoef = norm(D1(:, i));
        D1(:, i) = D1(:, i) / normCoef;
        D2(:, i) = D2(:, i) / normCoef;
    end
end
