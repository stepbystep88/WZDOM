function newResults = bsReOrganizeInvResults(invResults)
    nResults = length(invResults);
    
    types = {'seismic', 'ip', 'vp', 'vs', {'rho', 'density'}, {'vp_vs', 'vpvs_ratio'}, {'possion'}};
    newResults = {};
    for i = 1 : length(types)
        
        type = types{i};
        
        newProfile = [];
        
        for j = 1 : nResults
            profile = invResults{j};
            

            [index, rtype, rdata, rname] = bsGetIndexOfStrInCell(type, ...
                profile.type, profile.data, profile.name);
            
            if index > 0
                if isempty(newProfile)
                    newProfile = profile;
                    newProfile.type = {};
                    newProfile.data = {};
                    newProfile.name = {};
                end
            
                newProfile.type = [newProfile.type, rtype];
                newProfile.data = [newProfile.data, rdata];
                newProfile.name = [newProfile.name, rname];
            end
        end
        
        if ~isempty(newProfile)
            newResults = [newResults, newProfile];
        end
            
    end
    
end

function [index, rtype, rdata, rname] = bsGetIndexOfStrInCell(type, types, data, name)
    index = -1;
    rtype = [];
    rdata = [];
    rname = [];
    
    types = lower(types);
    type = lower(type);
    
    if iscell(types)
        for i = 1 : length(types)
            if bsStrCmp(type, types{i})
                index = i;
                rtype = types{i};
                rdata = data{i};
                
                if iscell(name)
                    rname = name{i};
                else
                    rname = name;
                end
                
                break;
            end
        end
    else
        if bsStrCmp(type, types)
            index = 1;
            rtype = types;
            rname = name;
            rdata = data;
        end
    end
end

function cmp = bsStrCmp(s1, s2)
    cmp = false;
    if iscell(s1)
        for i = 1 : length(s1)
            if strcmp(s1{i}, s2)
                cmp = true;
            end
        end
    else
        cmp = strcmp(s1, s2);
    end
end