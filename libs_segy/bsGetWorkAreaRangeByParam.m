function [rangeInline, rangeCrossline] = bsGetWorkAreaRangeByParam(GInvParam)
%% read the range of inline and crossline of a segy file
%
% Programmed by: Bin She (Email: bin.stepbystep@gmail.com)
% Programming dates: Nov 2019
% -------------------------------------------------------------------------
%     

    if isfield(GInvParam, 'postSeisData') ...
        && exist(GInvParam.postSeisData.segyFileName, 'file')
        
        [rangeInline, rangeCrossline] ...
            = bsGetWorkAreaRange(GInvParam.postSeisData.segyInfo, ...
                GInvParam.postSeisData.segyFileName);
    else
        [rangeInline, rangeCrossline] ...
            = bsGetWorkAreaRange(GInvParam.preSeisData.segyInfo, ...
                GInvParam.preSeisData.segyFileName);
    end                                        
end