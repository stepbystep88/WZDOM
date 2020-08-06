function [outLogs] = bsGetPairOfInvAndWell(GInvParam, timeLine, wellLogs, invResults, dataIndex, options)
    nWell = length(wellLogs);
    
    wells = cell2mat(wellLogs);
    wellInIds = [wells.inline];
    wellCrossIds = [wells.crossline];
    
    usedTimeLine = timeLine{GInvParam.usedTimeLineId};
    wellHorizonTimes = bsGetHorizonTime(usedTimeLine, ...
        wellInIds, wellCrossIds);
    indexInWellData = GInvParam.indexInWellData;
        GInvParam.upNum = GInvParam.upNum;
        GInvParam.downNum = GInvParam.downNum;
    outLogs = cell(1, nWell);
    startTime = wellHorizonTimes - GInvParam.upNum * GInvParam.dt;
    postSeisData = bsGetPostSeisData(GInvParam, wellInIds, wellCrossIds, startTime, GInvParam.upNum+GInvParam.downNum);
    
    % re-organize the data, and add the error data at the last column;
    for i = 1 : nWell
        
        wellInfo = wellLogs{i};
        

        subData = bsExtractWellDataByHorizon(...
                wellInfo.wellLog, ...
                wellHorizonTimes(i), ...
                [dataIndex, indexInWellData.time], ...
                indexInWellData.time, ...
                GInvParam.upNum, ...
                GInvParam.downNum, ...
                1);
       
        wellData = subData(:, 1);
        timeData = subData(:, 2);
        trueData = wellData;
        
        switch options.mode
            case 'full_freq'
                wellData = bsButtLowPassFilter(wellData, options.highCut);
                wellInfo.wellLog = [invResults(:, i), wellData];
            case 'low_high'
                % ´øÍ¨ÂË²¨
                if options.highCut < 1 && options.lowCut > 0
                    wellData = bsButtBandPassFilter(wellData, options.lowCut, options.highCut);
                end
                wellInfo.wellLog = [invResults(:, i), wellData];
                % µÍÍ¨ÂË²¨
%                 invResults(:, i) = bsButtLowPassFilter(invResults(:, i), options.lowCut);
            case 'seismic_high'
                wellData = bsButtBandPassFilter(wellData, options.lowCut, options.highCut);
                wellInfo.wellLog = [postSeisData(:, i), wellData];
            case 'residual'
                wellInfo.wellLog = [invResults(:, i), wellData - invResults(:, i)];
        end
        
        wellInfo.wellLog= [wellInfo.wellLog, trueData, timeData];
        outLogs{i} = wellInfo;
        
    end
    
end