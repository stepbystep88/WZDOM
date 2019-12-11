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
            profiles = {};
            for k = 1 : nResults
                if iscell(invResults{k}.data) && length(invResults{k}.data) > 1
                    
                    if length(invResults{k}.data) >= i
                        profile.data = invResults{k}.data{i};
                        profile.type = invResults{k}.type{i};

                        if iscell(invResults{k}.name)
                            profile.name = invResults{k}.name{i};
                        else
                            profile.name = invResults{k}.name;
                        end
                        
                        profiles = [profiles, profile];
                    end
                else
                    profile = bsRemoveCell(invResults{k});
                    profiles = [profiles, profile];
                end
            end
            
            if ~isempty(profiles)
                bsShowInvProfiles(GPreInvParam, GShowProfileParam, profiles, wellLogs, timeLine);
                set(gcf, 'position', [139          54        1690         797]);
            end
        end
        
    else
        bsShowInvProfiles(GPreInvParam, GShowProfileParam, profiles, wellLogs, timeLine);
    end
end
function profile = bsRemoveCell(profile)
    if iscell(profile.name)
        profile.name = profile.name{1};
    end
    if iscell(profile.data)
        profile.data = profile.data{1};
    end
    if iscell(profile.type)
        profile.type = profile.type{1};
    end
end