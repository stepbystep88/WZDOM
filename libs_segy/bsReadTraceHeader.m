function [traceHeader] = bsReadTraceHeader(GSegyInfo)
%% read the header information of one trace 
%
% Programmed by: Bin She (Email: bin.stepbystep@gmail.com)
% Programming dates: Nov 2019
% -------------------------------------------------------------------------
% Input
% GSegyInfo     basical infomation of the segy file
% 
% Output
% traceHeader   trace header information
% -------------------------------------------------------------------------
    
    bytes = fread(GSegyInfo.fid, 60, 'int32');                       
    
    traceHeader.inId = bytes(bsGetIntId(GSegyInfo.inlineId));
    traceHeader.crossId = bytes(bsGetIntId(GSegyInfo.crosslineId));
    traceHeader.X = bytes(bsGetIntId(GSegyInfo.xCoordId));
    traceHeader.Y = bytes(bsGetIntId(GSegyInfo.yCoordId));
    traceHeader.offset = bytes(bsGetIntId(GSegyInfo.offsetId));
    
    traceHeader.fullInfo = bytes;
end