function [GInvParam, wellLogs] = bsPreGetOtherAttributesOfWelllogs(GInvParam, wellLogs)
    for i = 1 : length(wellLogs)
        wellData = wellLogs{i}.wellLog;
        vpIndex = GInvParam.indexInWellData.vp;
        vsIndex = GInvParam.indexInWellData.vs;
        vp = wellData(:, vpIndex);
        vs = wellData(:, vsIndex);
        rho = wellData(:, GInvParam.indexInWellData.rho);
        
        vp_vs = bsGetVp_Vs(vp, vs);
        possion = bsGetPossion(vp, vs);
        brittleness = bsGetBrittleness(vp, vs, rho);
        
        toc = bsGetTOC(vp, vs, rho);
        
        
        wellLogs{i}.wellLog = [wellData, vp_vs, possion, brittleness, toc];
    end

    nAtt = size(wellData, 2);
    GInvParam.indexInWellData.vpvs_ratio = nAtt + 1;
    GInvParam.indexInWellData.possion = nAtt + 2;
    GInvParam.indexInWellData.brittleness = nAtt + 3;
    GInvParam.indexInWellData.toc = nAtt + 4;
    
end