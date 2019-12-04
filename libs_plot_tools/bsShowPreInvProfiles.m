function bsShowPreInvProfiles(GPreInvParam, GShowProfileParam, invResults, wellLogs, timeLine)
    nResults = length(invResults);
    
    profiles = invResults;
    hasMultiAtts = false;
    
    maxNumAtt = -1;
    for i = 1 : nResults
        profile = invResults{i};
        
        if iscell(profile.data)
            hasMultiAtts = true;
            if length(profile.data) > maxNumAtt
                maxNumAtt = length(profile.data);
            end
        end
    end
    
    if hasMultiAtts
        for i = 1 : maxNumAtt
            for k = 1 : nResults
                if iscell(invResults{k}.data)
                    profiles{k}.data = invResults{k}.data{i};
                    profiles{k}.type = invResults{k}.type{i};
                end
            end
            bsShowInvProfiles(GPreInvParam, GShowProfileParam, profiles, wellLogs, timeLine);
        end
        
    else
        bsShowInvProfiles(GPreInvParam, GShowProfileParam, profiles, wellLogs, timeLine);
    end
end