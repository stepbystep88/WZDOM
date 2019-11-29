function bsShowPreInvProfiles(GPreInvParam, GShowProfileParam, invResults, wellLogs)
    nResults = length(invResults);
    
    profiles = invResults;
    hasMultiAtts = false;
    
    for i = 1 : nResults
        profile = invResults{i};
        
        if iscell(profile.data)
            hasMultiAtts = true;
        end
    end
    
    if hasMultiAtts
        for i = 1 : 3
            for k = 1 : nResults
                if iscell(invResults{k}.data)
                    profiles{k}.data = invResults{k}.data{i};
                    profiles{k}.type = invResults{k}.type{i};
                end
            end
            bsShowInvProfiles(GPreInvParam, GShowProfileParam, profiles, wellLogs);
        end
        
    else
        bsShowInvProfiles(GPreInvParam, GShowProfileParam, profiles, wellLogs);
    end
end