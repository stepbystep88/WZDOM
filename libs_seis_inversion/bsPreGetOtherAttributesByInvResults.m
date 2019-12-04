function invResults = bsPreGetOtherAttributesByInvResults(invResults)
    
    vp = [];
    vs = [];
    rho = [];
    
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
    
end

function vp_vs = bsGetVp_Vs(vp, vs)
    vp_vs = vp ./ vs;
end

function possion = bsGetPossion(vp, vs)
    vp_vs = (vp ./ vs).^2;
    possion = (0.5*vp_vs - 1)./(vp_vs - 1);
end

function str = bsJointStr(strs, dchar)
    str = '';
    for i = 1 : length(strs)
        str = [str, dchar, strs{i}];
    end
end