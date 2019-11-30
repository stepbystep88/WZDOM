function [rangeInline, rangeCrossline] = bsGetWorkAreaRange(GSegyInfo, fileName)
%% read the range of inline and crossline of a segy file
%
% Programmed by: Bin She (Email: bin.stepbystep@gmail.com)
% Programming dates: Nov 2019
% -------------------------------------------------------------------------
%     

    GSegyInfo = bsReadVolHeader(fileName, GSegyInfo);        
    
    rangeInline = zeros(1, 2);
    rangeCrossline = zeros(1, 2);
    
    trHeader = bsReadTraceHeader(GSegyInfo);
    rangeCrossline(1) = trHeader.crossId;
    rangeInline(1) = trHeader.inId;
    
    volHeader = GSegyInfo.volHeader;
    fseek(GSegyInfo.fid, 3600 + (volHeader.traceNum-1)*(240+volHeader.sizeTrace), -1);
    trHeader = bsReadTraceHeader(GSegyInfo);
    rangeCrossline(2) = trHeader.crossId;
    rangeInline(2) = trHeader.inId;
    
    fclose(GSegyInfo.fid);                                        
end