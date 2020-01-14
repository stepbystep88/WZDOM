function [trData, sampNum, startTime] = bsReadTracesByIdsAndHorizons(fileName, GSegyInfo, inIds, crossIds, ...
    upHorizon, downHorizon, dt, valueRange)
%% read traces from a segy file with given inline and crossline ids
%
% Programmed by: Bin She (Email: bin.stepbystep@gmail.com)
% Programming dates: Nov 2019
% -------------------------------------------------------------------------
    startTime = bsGetHorizonTime(upHorizon, inIds, crossIds, 1);
    endTime = bsGetHorizonTime(downHorizon, inIds, crossIds, 1);
    
    sampNum = max( ceil((endTime - startTime) / dt) );
    
    trData = bsReadTraces(fileName, GSegyInfo, inIds, crossIds, ...
        startTime, sampNum, dt, valueRange);
end

function [trData, GSegyInfo] = bsReadTraces(fileName, GSegyInfo, inIds, crossIds, ...
    startTime, sampNum, dt, valueRange)
            
    GSegyInfo = bsReadVolHeader(fileName, GSegyInfo);        
    sizeTrace = GSegyInfo.volHeader.sizeTrace+240;    

    [~, trNum] = size(inIds);                            
    
    trData = zeros(sampNum, trNum);
    % the start position of a trace that we extract the data from
    startPos = bsGetT0Pos(GSegyInfo, startTime, dt);

    
    for i = 1 : trNum
        if mod(i, 10000) == 0
            % print information
            fprintf('Reading %d%% data from segy file %s...\n', round(i/trNum*100), fileName);
        end
        
        if i == 1
            index = bsIndexOfTraceSetOnInIdAndCrossId(GSegyInfo, inIds(i), crossIds(i));
        else
            index = bsIndexOfTraceSetOnInIdAndCrossId(GSegyInfo, inIds(i), crossIds(i), index);
        end
        
        if index > 0
            fseek(GSegyInfo.fid, 3600+sizeTrace*(index-1), -1);
            trHeader = bsReadTraceHeader(GSegyInfo);
            index = index + 1;
        else
            warning('Trace inline=%d, crossline=%d can not be found in file %s', inIds(i), crossIds(i), fileName);
            continue;
        end
        data = bsReadTraceData(GSegyInfo);
        
        iPos = startPos(i);
        tmp = data(iPos+1 : iPos+sampNum);
        
        % get the index of valid data
        valid = find(tmp >= valueRange(1) & tmp <= valueRange(2));
        
        if valid(1) < 1
            tmp(1:valid(1)-1) = tmp(valid(1));
        end
        
        if valid(end) < sampNum
            tmp(valid(end):end) = tmp(valid(end));
        end
        
        trData(:, i) = tmp;
        
    end

    fclose(GSegyInfo.fid);    
    

end