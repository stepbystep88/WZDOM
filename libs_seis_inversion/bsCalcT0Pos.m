function pos = bsCalcT0Pos(GInvParam, GSegyInfo, time)
    if( GSegyInfo.isPosZero )
        pos = ones(size(time)) * 0;
    else
        pos = round((time - GSegyInfo.t0) / GInvParam.dt) - GInvParam.upNum;
    end
end