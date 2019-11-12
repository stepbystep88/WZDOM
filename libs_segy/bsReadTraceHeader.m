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
    
    traceHeader.fullInfo = fread(GSegyInfo.fid, 60, 'int32');                       
    
    traceHeader.inId = traceHeader.fullInfo(GSegyInfo.inlineId, 1);                 % inline id
    traceHeader.crossId = traceHeader.fullInfo(GSegyInfo.crosslineId, 1);           % crossline id
    traceHeader.X = traceHeader.fullInfo(GSegyInfo.xCoordId, 1);
    traceHeader.Y = traceHeader.fullInfo(GSegyInfo.yCoordId, 1);
    traceHeader.offset = traceHeader.fullInfo(GSegyInfo.offsetId, 1);               % offsetId
end