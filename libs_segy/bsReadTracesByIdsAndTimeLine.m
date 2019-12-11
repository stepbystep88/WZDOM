function [trData, GSegyInfo] = bsReadTracesByIdsAndTimeLine(fileName, GSegyInfo, inIds, crossIds, ...
    usedTimeLine, upNum, downNum, dt)
%% read traces from a segy file with given inline and crossline ids
%
% Programmed by: Bin She (Email: bin.stepbystep@gmail.com)
% Programming dates: Nov 2019
% -------------------------------------------------------------------------
    horizonTime = bsGetHorizonTime(usedTimeLine, inIds, crossIds, 1);
    startTime = horizonTime - upNum * dt;
    sampNum = upNum + downNum;
    
    trData = bsReadTracesByIds(fileName, GSegyInfo, inIds, crossIds, startTime, sampNum, dt);
end