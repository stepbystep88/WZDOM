function [wellLogs, timeLine, GInvParam] = bsPreGenWelllogs(GInvParam, profileData, well_ids, t0)
    nWell = length(well_ids);
    
    wellLogs = cell(1, nWell);
    [sampNum, traceNum] = size(profileData);
    
    time = (0:sampNum-1)'*GInvParam.dt + t0;
    
    GInvParam.indexInWellData.ip = 6;
    GInvParam.indexInWellData.vp = 2;
    GInvParam.indexInWellData.vs = 3;
    GInvParam.indexInWellData.rho = 4;
    GInvParam.indexInWellData.depth = 1;
    GInvParam.indexInWellData.time = 5;
        
    for i = 1 : nWell
        t.name = sprintf('#%d', i);
        t.inline = 1;
        t.crossline = well_ids(:, i);
        
        k = well_ids(:, i);
        
        vp = profileData(:, k, 2);
        vs = profileData(:, k, 3);
        rho = profileData(:, k, 4);
        ip = vp .* rho;
        
        depth = bsGetDepth(vp, GInvParam.dt);
        
        
        t.wellLog = [depth, vp, vs, rho, time, ip];
        
        wellLogs{i} = t;
    end
    
    timeLine = {ones(1, traceNum) * t0};
end