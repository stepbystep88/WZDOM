function ids = bsGetIdsFromWelllogs(welllogs, names)
    
    ids = zeros(1, length(names));
    
    for k = 1 : length(names)
        
        name = names{k};
        
        for i = 1 : length(welllogs)
            
            wellInfo = welllogs{i};
            if isfield(wellInfo, 'name') && strcmp(wellInfo.name, name)
                ids(k) = i;
                break;
            elseif isfield(wellInfo, 'wellName') && strcmp(wellInfo.wellName, name)
                ids(k) = i;
            end
            
        end
    end
    
    ids(find(ids==0)) = [];
end