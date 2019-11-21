function bsWriteTrace(GSegyInfo, traceHeader, data)
%% write one trace into segy file
%
% Programmed by: Bin She (Email: bin.stepbystep@gmail.com)
% Programming dates: Nov 2019
% -------------------------------------------------------------------------
% Input
% GSegyInfo     basical infomation of the segy file
% traceHeader   trace header information
% data          trace data to save
% -------------------------------------------------------------------------

    if(~GSegyInfo.isSaveCoordInfo)
        traceHeader.fullInfo(bsGetIntId(GSegyInfo.xCoordId)) = 0;
        traceHeader.fullInfo(bsGetIntId(GSegyInfo.yCoordId)) = 0;
    end
    
    traceHeader.fullInfo(bsGetIntId(GSegyInfo.inlineId)) = traceHeader.inId;                         
    traceHeader.fullInfo(bsGetIntId(GSegyInfo.crosslineId)) = traceHeader.crossId;                   
    traceHeader.fullInfo(bsGetIntId(GSegyInfo.offsetId)) = traceHeader.offset;                       
    traceHeader.fullInfo(bsGetIntId(GSegyInfo.traceSampNumId)) = length(data);
    
    fwrite(GSegyInfo.fid, traceHeader.fullInfo, 'int32');
    
    fwrite(GSegyInfo.fid, data, 'float32');
end

