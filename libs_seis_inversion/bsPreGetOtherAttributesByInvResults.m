function [invResults, GInvParam, wellLogs] = bsPreGetOtherAttributesByInvResults(invResults, GInvParam, wellLogs)
    
    vp = [];
    vs = [];
    
    for i = 1 : length(invResults)
        profile = invResults{i};
        nData = length(profile.data);
        
        if iscell(profile.data)
            for j = 1 : nData
                switch lower(profile.type{j})
                    case 'vp'
                        vp = profile.data{j};
                    case 'vs'
                        vs = profile.data{j};
                    case {'rho', 'density'}
                        rho = profile.data{j};
                end
            end
            
            vp_vs = bsGetVp_Vs(vp, vs);
            possion = bsGetPossion(vp, vs);

            profile.data{nData+1} = vp_vs;
            profile.type{nData+1} = 'vp_vs';
            profile.data{nData+2} = possion;
            profile.type{nData+2} = 'possion';
        end
        invResults{i} = profile;
    end
    
    if nargin > 1
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
    else
        GInvParam = [];
        wellLogs = [];
    end
    
end



