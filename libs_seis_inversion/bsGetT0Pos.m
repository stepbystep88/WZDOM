function pos = bsGetT0Pos(GSegyInfo, startTime, dt)
    if( GSegyInfo.isPosZero )
        pos = ones(size(startTime)) * 0;
    else
        pos = round((startTime - GSegyInfo.t0) / dt);
    end
end