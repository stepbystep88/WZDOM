function [outputData, highData] = ...
    bsPostReBuildInterpolation(GInvParam, wellLogs, inputData, inIds, crossIds, options)

    [sampNum, traceNum] = size(inputData);
    
        % tackle the inverse task
    outputData = zeros(sampNum, traceNum);
    highData = zeros(sampNum, traceNum);
    
    dt = GInvParam.dt;
    
    wells = cell2mat(wellLogs);
    wellInIds = [wells.inline];
    wellCrossIds = [wells.crossline];

    wellData = [];
    for i = 1 : length(wellLogs)
        wellData = [wellData, wellLogs{i}.wellLog(:, 2)];
    end
    
    op.p = 2;
    op.nPointsUsed = 4;
    op.filtCoef = 1;
    
    [weights, indexies] = bsGetWeightByIDW(inIds, crossIds, wellInIds, wellCrossIds, op);
    % 给未高分辨率的道插值
    highData = bsInterpolate3DData(length(inIds), wellData, weights, indexies);
    
    % parallel computing
    for iTrace = 1 : traceNum
        outputData(:, iTrace) = bsHandleOneTrace(inputData(:, iTrace), highData(:, iTrace), options, dt);
    end

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

  