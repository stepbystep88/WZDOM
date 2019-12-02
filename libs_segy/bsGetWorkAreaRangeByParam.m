function [rangeInline, rangeCrossline, fileName, segyInfo] ...
    = bsGetWorkAreaRangeByParam(GInvParam)
%% read the range of inline and crossline of a segy file
%
% Programmed by: Bin She (Email: bin.stepbystep@gmail.com)
% Programming dates: Nov 2019
% -------------------------------------------------------------------------
    if isfield(GInvParam, 'postSeisData') ...
        && exist(GInvParam.postSeisData.segyFileName, 'file')
        
        fileName = GInvParam.postSeisData.segyFileName;
        segyInfo = GInvParam.postSeisData.segyInfo;
        
    else
        fileName = GInvParam.preSeisData.segyFileName;
        segyInfo = GInvParam.preSeisData.segyInfo;
    end    
    
    [rangeInline, rangeCrossline] = bsGetWorkAreaRange(segyInfo, fileName);
end