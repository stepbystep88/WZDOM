function [inIds, crossIds] = bsBuildInitModel(GInvParam, timeLine, wellLogs, varargin)
%% Build initial model and save the result as segy file
% Programmed by: Bin She (Email: bin.stepbystep@gmail.com)
% Programming dates: Dec 2019
% -------------------------------------------------------------------------

    p = inputParser;
    
    [rangeInline, rangeCrossline, fileName, segyInfo] ...
        = bsGetWorkAreaRangeByParam(GInvParam);
    
    addParameter(p, 'filtCoef', 0.1);
    addParameter(p, 'title', bsGetCurentDateAsString());
    addParameter(p, 'prestack', 1);
    addParameter(p, 'dstPath', sprintf('%s/model/', GInvParam.modelSavePath));
    addParameter(p, 'rangeInline', rangeInline);
    addParameter(p, 'rangeCrossline', rangeCrossline);
    addParameter(p, 'nMostUseWells', 4);
    addParameter(p, 'p', -2);
    
    p.parse(varargin{:});  
    options = p.Results;
    
    [inIds, crossIds] = bsGetCDPsByRange(options.rangeInline, options.rangeCrossline);
    
    wells = cell2mat(wellLogs);
    wellInIds = [wells.inline];
    wellCrossIds = [wells.crossline];
    
    usedTimeLine = timeLine{GInvParam.usedTimeLineId};
    wellHorizonTimes = bsCalcHorizonTime(usedTimeLine, ...
        wellInIds, wellCrossIds);
    horizonTimes = bsCalcHorizonTime(usedTimeLine, inIds, crossIds, ...
            GInvParam.isParallel, GInvParam.numWorkers);
        
    assert(length(inIds) == length(crossIds),...
        'The length of inIds must be the same as that of crossIds.');
    
    if options.prestack 
        dataIndex = [...
            GInvParam.indexInWellData.vp, ...
            GInvParam.indexInWellData.vs, ...
            GInvParam.indexInWellData.rho];
        type = {'vp', 'vs', 'density'};
    else
        dataIndex = GInvParam.indexInWellData.Ip;
        type = {'ip'};
    end
    
    [weights, indexies] = bsGetWeightByIDW(inIds, crossIds, ...
        wellInIds, wellCrossIds, options.nMostUseWells, -options.p);
    nTrace = length(inIds);
    res.inIds = inIds;
    res.crossIds = crossIds;
    res.horizon = horizonTimes;
    res.upNum = GInvParam.upNum;
    res.dt = GInvParam.dt;
    
    for i = 1 : length(dataIndex)
        wellData = bsGetWellData(GInvParam, wellLogs, wellHorizonTimes, dataIndex(i), options.filtCoef);
        
        fprintf('Interpolating the %s data by calculated weights...\n', type{i});
        data = bsInterpolate3DData(nTrace, wellData, weights, indexies);
        
        dstFileName = sprintf('%s/%s_%s_flitCoef_%.2f.sgy', ...
            options.dstPath, type{i}, options.title, options.filtCoef);
        
        bsWriteInvResultIntoSegyFile(res, data, fileName, segyInfo, dstFileName);
    end
    
end

function wellData = bsGetWellData(GInvParam, wellLogs, wellHorizonTimes, dataIndex, filtCoef)

    sampNum = GInvParam.upNum + GInvParam.downNum;
    wellNum = length(wellLogs);
    
    wellData = zeros(sampNum, wellNum);
    
    for i = 1 : wellNum
        
        wellData(:, i) = bsExtractWellDataByHorizon(...
            wellLogs{i}.wellLog, ...
            wellHorizonTimes(i), ...
            dataIndex, ...
            GInvParam.indexInWellData.time, ...
            GInvParam.upNum, ...
            GInvParam.downNum, ...
            filtCoef);
        
    end
end

function [weights, indexies] = bsGetWeightByIDW(inIds, crossIds, ...
    wellInIds, wellCrossIds, nMostUseWells, e)

    nTrace = length(inIds);
    
    fprintf('Calculating the weight information...\n');
    
     % at most use 4 wells for interpolation
    nUseWell = min(nMostUseWells, length(wellInIds));
    indexies = zeros(nTrace, nUseWell);
    weights = zeros(nUseWell, nTrace);
    
    for i = 1 : nTrace
        D = sqrt((inIds(i)-wellInIds).^2 + (crossIds(i)-wellCrossIds).^2);
        [~, index] = bsMinK(D, nUseWell); 
        indexies(i, :) = index';
        
        [minD, minIndex] = min(D);
        if minD == 0
            weights(1, i) = 1;
        else
            weights(:, i) = D(index).^e;
        end
        weights(:, i) = weights(:, i) / sum(weights(:, i));
        
        if mod(i, 10000) == 0
            fprintf('Calculating the weight information of trace %d/%d...\n', i, nTrace);
        end
    end
    
end

function data = bsInterpolate3DData(nTrace, wellData, weights, indexies)

    sampNum = size(wellData, 1);
    data = zeros(sampNum, nTrace);
    
    for i = 1 : nTrace
        weight = weights(:, i);
        usedWellData = wellData(:, indexies(i, :));
        
        V = usedWellData * weight;
        if i == 442
            i;
        end
        
        data(:, i) = V;
        
        if mod(i, 10000) == 0
            fprintf('Interpolating the model of trace %d/%d...\n', i, nTrace);
        end
    end
end
