function bsWriteInvResultIntoSegyFile(res, sourceFileName, sourceSegyInfo, dstFileName, upNum, dt)
    
    data = res.data;
    inIds = res.inIds;
    crossIds = res.crossIds;
    horizon = res.horizon;
    
   
    [newProfileData, minTime] = bsHorizonRestoreData(data, horizon, upNum, dt, sourceSegyInfo.t0, 1);
    newProfileData(isnan(newProfileData)) = -99999;
    
    bsWriteTracesByRefFileAndIds(sourceFileName, dstFileName, sourceSegyInfo, newProfileData, inIds, crossIds);
    
%     outSegInfo = bsWriteTracesByIds(dstFileName, sourceSegyInfo, newProfileData, inIds, crossIds, traceHeader);
    
    fprintf('\nStart time of file %s is %.2f \n\n', dstFileName, minTime);
end