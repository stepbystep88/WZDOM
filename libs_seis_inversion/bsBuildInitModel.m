function [inIds, crossIds, GInvParam] = bsBuildInitModel(GInvParam, timeLine, wellLogs, varargin)
%% Build initial model and save the result as segy file
% Programmed by: Bin She (Email: bin.stepbystep@gmail.com)
% Programming dates: Dec 2019
% -------------------------------------------------------------------------

    p = inputParser;
    
    [rangeInline, rangeCrossline, fileName, segyInfo] ...
        = bsGetWorkAreaRangeByParam(GInvParam);
    
    addParameter(p, 'filtCoef', 0.1);
    addParameter(p, 'title', '');
    addParameter(p, 'prestack', 1);
    addParameter(p, 'dstPath', sprintf('%s/model/', GInvParam.modelSavePath));
    addParameter(p, 'rangeInline', rangeInline);
    addParameter(p, 'rangeCrossline', rangeCrossline);
    addParameter(p, 'nPointsUsed', 4);
    addParameter(p, 'p', 1.2);
    addParameter(p, 'expandNum', 30);
    
    p.parse(varargin{:});  
    options = p.Results;
    
    [inIds, crossIds] = bsGetCDPsByRange(options.rangeInline, options.rangeCrossline);
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
    
    [isExist, GInvParam] = checkExists(GInvParam, type, dataIndex, options, segyInfo);
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
    
    for i = 1 : length(dataIndex)
        wellData = bsGetWellData(GInvParam, wellLogs, wellHorizonTimes, dataIndex(i), options.filtCoef);
        
        fprintf('Interpolating the %s data by calculated weights...\n', type{i});
        data = bsInterpolate3DData(nTrace, wellData, weights, indexies);
        
        dstFileName = bsGetDstFileName(type{i}, options);
        bsWriteInvResultIntoSegyFile(res, data, fileName, segyInfo, dstFileName);
        
    end
    
    GInvParam.upNum = GInvParam.upNum - options.expandNum;
    GInvParam.downNum = GInvParam.downNum - options.expandNum;
    
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

function [isExist, GInvParam] = checkExists(GInvParam, type, dataIndex, options, segyInfo)
    isExist = true;
    for j = 1 : length(dataIndex)

        dstFileName = bsGetDstFileName(type{j}, options);

        if exist(dstFileName, 'file')
            warning('Initial model %s exists already!', dstFileName);
        else
            isExist = false;
        end

        % update GInvParam;
        switch lower(type{j})
            case 'vp'
                GInvParam.initModel.vp.segyFileName = dstFileName;
                GInvParam.initModel.vp.segyInfo = segyInfo;
            case 'vs'
                GInvParam.initModel.vs.segyFileName = dstFileName;
                GInvParam.initModel.vs.segyInfo = segyInfo;
            case 'density'  
                GInvParam.initModel.rho.segyFileName = dstFileName;
                GInvParam.initModel.rho.segyInfo = segyInfo;
            case 'ip'
                GInvParam.initModel.segyFileName = dstFileName;
                GInvParam.initModel.segyInfo = segyInfo;
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
