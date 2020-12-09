function [invResults, GInvParam, wellLogs] = bsPreGetOtherAttributesByInvResults(invResults, GInvParam, wellLogs, varargin)
    
    p = inputParser;
    addParameter(p, 'atts', {'vp_vs', 'possion', 'brittleness', 'toc'});
    p.parse(varargin{:});  
    options = p.Results;
    
    vp = [];
    vs = [];
    
    if length(invResults) < 2
        return;
    end
    
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
            
            for k = 1 : length(options.atts)
                switch lower(options.atts{k})
                    case 'vp_vs'
                        data = bsGetVp_Vs(vp, vs);
                        type = 'vp_vs';
                    case 'possion'
                        data = bsGetPossion(vp, vs);
                        type = 'possion';
                    case {'brittleness', 'cuixing', 'cui_xing'}
                        data = bsGetBrittleness(vp, vs, rho);
                        type = 'brittleness';
                    case 'toc'
                        data = bsGetTOC(vp, vs, rho);
                        type = 'toc';
                end
                
                profile.data{nData+k} = data;
                profile.type{nData+k} = type;
            end
            
        end
        invResults{i} = profile;
    end
    
    if nargin > 1
         [GInvParam, wellLogs] = bsPreGetOtherAttributesOfWelllogs(GInvParam, wellLogs);
    else
        GInvParam = [];
        wellLogs = [];
    end
    
end



