function bsWriteInvResultIntoSegyFile(res, sourceFileName, sourceSegyInfo, dstFileName, upNum, dt)
    
    data = res.data;
    inIds = res.inIds;
    crossIds = res.crossIds;
    horizon = res.horizon;
    
    [sourceSegyInfo] = bsReadVolHeader(sourceFileName, sourceSegyInfo);
    [traceHeader] = bsReadTraceHeader(sourceSegyInfo);
    fclose(sourceSegyInfo.fid);
    
    data(isnan(data)) = -9999;
    [newProfileData, minTime] = bsHorizonRestoreData(data, horizon, upNum, dt, 1);
    outSegInfo = bsWriteTracesByIds(dstFileName, sourceSegyInfo, newProfileData, inIds, crossIds, traceHeader);
    
    fprintf('\nStart time of file %s is %.2f \n', dstFileName, minTime);
end