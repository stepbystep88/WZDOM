function [inIds, crossIds, GInvParam, dstFileNames, segyInfo] = bsBuildInitModel(GInvParam, timeLine, wellLogs, varargin)
%% Build initial model and save the result as segy file
% Programmed by: Bin She (Email: bin.stepbystep@gmail.com)
% Programming dates: Dec 2019
% -------------------------------------------------------------------------
    
    p = inputParser;
    
    [rangeInline, rangeCrossline, fileName, segyInfo] ...
        = bsGetWorkAreaRangeByParam(GInvParam);
    
    addParameter(p, 'filtCoef', 0.1);
    addParameter(p, 'title', '');
    addParameter(p, 'dstPath', sprintf('%s/model/', GInvParam.modelSavePath));
    addParameter(p, 'rangeInline', rangeInline);
    addParameter(p, 'rangeCrossline', rangeCrossline);
    addParameter(p, 'inIds', []);
    addParameter(p, 'crossIds', []);
    addParameter(p, 'nPointsUsed', 10);
    addParameter(p, 'p', 2);
    addParameter(p, 'expandNum', 30);
    addParameter(p, 'isRebuild', 0);
    addParameter(p, 'dataIndex', []);
    addParameter(p, 'type', []);
    
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
    
    assert(length(inIds) == length(crossIds),...
        'The length of inIds must be the same as that of crossIds.');
    
    if isempty(options.dataIndex)
        switch lower(GInvParam.flag)
        case {'prestack', 'pre-stack'}
            dataIndex = [...
                GInvParam.indexInWellData.vp, ...
                GInvParam.indexInWellData.vs, ...
                GInvParam.indexInWellData.rho];
            type = {'vp', 'vs', 'density'};
        case {'poststack', 'post-stack'}
            dataIndex = GInvParam.indexInWellData.ip;
            type = {'ip'};
        end
    else
        dataIndex = options.dataIndex;
        type = options.type;
    end
    
    [isExist, GInvParam, dstFileNames] = checkExists(GInvParam, type, dataIndex, options, segyInfo);
    if isExist
        return;
    end
    mkdir(options.dstPath);
    
    GInvParam.upNum = GInvParam.upNum + options.expandNum;
    GInvParam.downNum = GInvParam.downNum + options.expandNum;
    
    wells = cell2mat(wellLogs);
    wellInIds = [wells.inline];
    wellCrossIds = [wells.crossline];
    
    usedTimeLine = timeLine{GInvParam.usedTimeLineId};
    wellHorizonTimes = bsGetHorizonTime(usedTimeLine, ...
        wellInIds, wellCrossIds);
    horizonTimes = bsGetHorizonTime(usedTimeLine, inIds, crossIds, ...
            GInvParam.isParallel, GInvParam.numWorkers);
        
    % calculate the weight information
    [weights, indexies] = bsGetWeightByIDW(inIds, crossIds, ...
        wellInIds, wellCrossIds, options);
    nTrace = length(inIds);
    res.inIds = inIds;
    res.crossIds = crossIds;
    res.horizon = horizonTimes;
    res.upNum = GInvParam.upNum;
    res.dt = GInvParam.dt;
    
    nData = length(dataIndex);
    for i = 1 : nData
        wellData = bsGetWellData(GInvParam, wellLogs, wellHorizonTimes, dataIndex(i), options.filtCoef);
        
        fprintf('Interpolating the %s data by calculated weights...\n', type{i});
        data = bsInterpolate3DData(nTrace, wellData, weights, indexies);
        
        dstFileName = bsGetDstFileName(type{i}, options);
        bsWriteInvResultIntoSegyFile(res, data, fileName, segyInfo, dstFileName);
    end
    
    GInvParam.upNum = GInvParam.upNum - options.expandNum;
    GInvParam.downNum = GInvParam.downNum - options.expandNum;
    
end
   


function [isExist, GInvParam, dstFileNames] = checkExists(GInvParam, type, dataIndex, options, segyInfo)
    isExist = true;
    dstFileNames = cell(1, length(dataIndex));
    
    for j = 1 : length(dataIndex)

        dstFileName = bsGetDstFileName(type{j}, options);
        dstFileNames{j} = dstFileName;
        
        if exist(dstFileName, 'file') && ~options.isRebuild
            warning('Initial model %s exists already!', dstFileName);
        else
            isExist = false;
        end

        % update GInvParam;
        switch lower(type{j})
            case 'vp'
                GInvParam.initModel.vp.fileName = dstFileName;
                GInvParam.initModel.vp.segyInfo = segyInfo;
            case 'vs'
                GInvParam.initModel.vs.fileName = dstFileName;
                GInvParam.initModel.vs.segyInfo = segyInfo;
            case 'density'  
                GInvParam.initModel.rho.fileName = dstFileName;
                GInvParam.initModel.rho.segyInfo = segyInfo;
            case 'ip'
                GInvParam.initModel.ip.fileName = dstFileName;
                GInvParam.initModel.ip.segyInfo = segyInfo;
        end
    end
    GInvParam.initModel.mode = 'segy';
end
    
function dstFileName = bsGetDstFileName(type, options)
    dstFileName = sprintf('%s/%s_%s_flitCoef_%.2f_p_%.2f_inline_[%d_%d]_crossline_[%d_%d].sgy', ...
        options.dstPath, type, options.title, options.filtCoef, options.p, ...
        options.rangeInline(1), options.rangeInline(end), ...
        options.rangeCrossline(1), options.rangeCrossline(2));
end

function data = bsInterpolate3DData(nTrace, wellData, weights, indexies)

    sampNum = size(wellData, 1);
    data = zeros(sampNum, nTrace);
    
    for i = 1 : nTrace
        weight = weights(:, i);
        usedWellData = wellData(:, indexies(i, :));
        
        V = usedWellData * weight;
        
        data(:, i) = V;
        
        if mod(i, 10000) == 0
            fprintf('Interpolating the model of trace %d/%d...\n', i, nTrace);
        end
    end
end
