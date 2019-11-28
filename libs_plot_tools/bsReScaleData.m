function [newProfileData, newDt, newTraceIds, newHorizon] = bsReScaleData(scaleFactor, profileData, dt, traceIds, horizon)
    if exist('scaleFactor', 'var') && scaleFactor > 1
        scaleFactor = round(scaleFactor);
        [sampNum, traceNum] = size(profileData);
        
        newDt = dt / scaleFactor;
        newTraceIds = linspace(traceIds(1), traceIds(end), traceNum * 1);
        [X, Y] = meshgrid(traceIds, linspace(1, sampNum, sampNum));
        [Xq,Yq] = meshgrid(newTraceIds, linspace(1, sampNum, sampNum * scaleFactor));
        
        newProfileData = interp2(X, Y, profileData, Xq, Yq,'cubic');
        newHorizon = interp1(traceIds, horizon, newTraceIds, 'spline');
    else
        newProfileData = profileData;
        newDt = dt;
        newTraceIds = traceIds;
        newHorizon = horizon;
    end
end