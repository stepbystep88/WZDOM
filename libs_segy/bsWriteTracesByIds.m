function GSegyInfo = bsWriteTracesByIds(fileName, GSegyInfo, trDatas, inIds, crossIds, trHeaders)
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

    GSegyInfo.sampNum = size(trDatas, 1);
    GOutSegyInfo = bsWriteVolHeader(fileName, GSegyInfo);
 
    if ~iscell(trHeaders) || length(trHeaders) ~= size(trDatas, 2)
        error('Header must be a struct or a cell with the same length of trace data.');
    end
        
    for i = 1 : size(trDatas, 2)
        if iscell(trHeaders)
            bsWriteTrace(GOutSegyInfo, trHeaders{i}, trDatas{i});
        else
            trHeaders.inId = inIds(i);
            trHeaders.crossId = crossIds(i);
            
            bsWriteTrace(GOutSegyInfo, trHeaders, trDatas{i});
        end
    end
    
    fclose(GOutSegyInfo.fid);
    
end
