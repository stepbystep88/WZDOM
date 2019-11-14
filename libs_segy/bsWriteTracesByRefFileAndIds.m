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
    
    for i = 1 : length(inIds)
        index = bsIndexOfTraceSetOnInIdAndCrossId(sourceSegyInfo, inIds(i), crossIds(i));
        
        if index > 0
            fseek(sourceSegyInfo.fid, 3600+sizeTrace*(index-1), -1);
            traceHeader = bsReadTraceHeader(sourceSegyInfo);
            bsWriteTrace(outSegyInfo, traceHeader, trDatas(:, i));
        end
    end
    
    fclose(sourceSegyInfo.fid);
    fclose(outSegyInfo.fid);
end
