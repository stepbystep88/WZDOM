function pos = bsCalcT0Pos(GPostInvParam, GSegyInfo, time)
    if( GSegyInfo.isPosZero )
        pos = 0;
    else
        pos = round((time - GSegyInfo.t0) / GPostInvParam.dt) - GPostInvParam.upNum;
    end
end