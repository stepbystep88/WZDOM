function [firstInline, firstCrossline, endInline, endCrossline] = bsGetWorkAreaRange(GSegyInfo, fileName)
%% read the range of inline and crossline of a segy file
%
% Programmed by: Bin She (Email: bin.stepbystep@gmail.com)
% Programming dates: Nov 2019
% -------------------------------------------------------------------------
%     

    GSegyInfo = bsReadVolHeader(fileName, GSegyInfo);        
    

    trHeader = bsReadTraceHeader(GSegyInfo);
    firstCrossline = trHeader.crossId;
    firstInline = trHeader.inId;
    
    volHeader = GSegyInfo.volHeader;
    fseek(GSegyInfo.fid, 3600 + (volHeader.traceNum-1)*(240+volHeader.sizeTrace), -1);
    trHeader = bsReadTraceHeader(GSegyInfo);
    endCrossline = trHeader.crossId;
    endInline = trHeader.inId;
    
    fclose(GSegyInfo.fid);                                        
end