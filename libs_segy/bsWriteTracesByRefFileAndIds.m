function bsWriteTracesByRefFileAndIds(refFileName, dstfileName, GSegyInfo, trDatas, inIds, crossIds)
%% write multiple traces into a segy file
%
% Programmed by: Bin She (Email: bin.stepbystep@gmail.com)
% Programming dates: Nov 2019
% -------------------------------------------------------------------------
% Input
% fileName          name of the segy file to save the data
% GSegyInfo         basical information of the segy file
% inIds             inline ids
% crossIds          crossline ids
% trDatas           a matrix saving the traceData of each 
% -------------------------------------------------------------------------

    [sourceSegyInfo] = bsReadVolHeader(refFileName, GSegyInfo);
    sizeTrace = sourceSegyInfo.volHeader.sizeTrace+240;    
    
    outSegyInfo = sourceSegyInfo;
    outSegyInfo.volHeader.sampNum = size(trDatas, 1);
    outSegyInfo = bsWriteVolHeader(dstfileName, outSegyInfo);
    
    trNum = length(inIds);
    index = -1;
    
    for i = 1 : trNum
        if mod(i, 1000) == 0
            % print information
            fprintf('Writing %d%% data into segy file %s...\n', round(i/trNum*100), dstfileName);
        end
        
        if i > 1
            [index, traceHeader] = bsMoveToNextTrace(index, inIds(i), crossIds(i));
        end
        
        if i == 1 || index < 0
            index = bsIndexOfTraceSetOnInIdAndCrossId(sourceSegyInfo, inIds(i), crossIds(i));
            fseek(sourceSegyInfo.fid, 3600+sizeTrace*(index-1), -1);
            traceHeader = bsReadTraceHeader(sourceSegyInfo);
        end
        
        if index > 0
            bsWriteTrace(outSegyInfo, traceHeader, trDatas(:, i));
        else
            error('Trace inline=%d, crossline=%d can not be found in file %s', inIds(i), crossIds(i), ref);
        end
    end
    
    fclose(sourceSegyInfo.fid);
    fclose(outSegyInfo.fid);
    
    function [index, traceHeader] = bsMoveToNextTrace(index, nextInId, nextCrossId)
        fseek(sourceSegyInfo.fid, 3600+sizeTrace*index, -1);
        traceHeader = bsReadTraceHeader(sourceSegyInfo);

        if traceHeader.inId == nextInId && traceHeader.crossId == nextCrossId
            index = index + 1;
        else
            index = -1;
            traceHeader = [];
        end
    end

end




