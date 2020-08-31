function [GInvParam, wellLogs] = bsPreGetOtherAttributesOfWelllogs(GInvParam, wellLogs)
    for i = 1 : length(wellLogs)
        wellData = wellLogs{i}.wellLog;
        vpIndex = GInvParam.indexInWellData.vp;
        vsIndex = GInvParam.indexInWellData.vs;
        vp = wellData(:, vpIndex);
        vs = wellData(:, vsIndex);
        vp_vs = bsGetVp_Vs(vp, vs);
        possion = bsGetPossion(vp, vs);
        wellLogs{i}.wellLog = [wellData, vp_vs, possion];
    end

    nAtt = size(wellData, 2);
    GInvParam.indexInWellData.vpvs_ratio = nAtt + 1;
    GInvParam.indexInWellData.possion = nAtt + 2;

end