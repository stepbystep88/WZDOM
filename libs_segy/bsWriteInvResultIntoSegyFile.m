function bsWriteInvResultIntoSegyFile(res, data, sourceFileName, sourceSegyInfo, dstFileName, isSort)
    
    inIds = res.inIds;
    crossIds = res.crossIds;
    horizon = res.horizon;
    
    if nargin <= 5
        isSort = 1;
    end
    
    if isSort
        % resort the ids 
        [ids, index] = sortrows([inIds', crossIds'], [1, 2]);
        index = index';

        inIds = ids(:, 1);
        crossIds = ids(:, 2);
        data = data(:, index);
        horizon = horizon(index);
    end
    
    [newProfileData, minTime] = bsHorizonRestoreData(data, horizon, res.upNum, res.dt, sourceSegyInfo.t0, 1);
    newProfileData(isnan(newProfileData)) = -99999;
    
    bsWriteTracesByRefFileAndIds(sourceFileName, dstFileName, sourceSegyInfo, newProfileData, inIds, crossIds);
    
%     outSegInfo = bsWriteTracesByIds(dstFileName, sourceSegyInfo, newProfileData, inIds, crossIds, traceHeader);
    
    fprintf('\nStart time of the sgy file %s is %.2fms \n\n', dstFileName, minTime);
end