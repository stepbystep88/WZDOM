function bsCheckIsValidWellData(GInvParam, wellLogs, varargin)
    p = inputParser;
    
    addParameter(p, 'range', struct(...
        'ip', [1000 20000], ...
        'vp', [1000 8000], ...
        'vs', [0 7000], ...
        'rho', [0 5], ...
        'vpvs_ratio', [0 5], ...
        'possion', [0 100] ...
    )); 
    p.parse(varargin{:});  
    options = p.Results;
    
    indexInWell = GInvParam.indexInWellData;
    fields = fieldnames(indexInWell); % cell
    
    for i = 1 : length(fields)
        name_i = fields{i};
        index = getfield(indexInWell, name_i);
        
        if index <= 0 || index > size(wellLogs{1}.wellLog, 2)
            continue;
        end
        
        for k = 1 : length(wellLogs)
            maxVal = max(wellLogs{k}.wellLog(:, index));
            minVal = min(wellLogs{k}.wellLog(:, index));
            
            if isempty(maxVal) || isempty(minVal)
                warning('The %s welllog data of %d-th well is empty', name_i, k);
                continue;
            end
            
            switch  name_i
                case 'vp' 
                    if maxVal > options.range.vp(2) || minVal < options.range.vp(1)
                        warning('The vp of %d-th well is not valid, please check its data range.', k);
                    end
                case 'vs' 
                    if maxVal > options.range.vs(2) || minVal < options.range.vs(1)
                        warning('The vs of %d-th well is not valid, please check its data range.', k);
                    end
                case 'ip' 
                    if maxVal > options.range.ip(2) || minVal < options.range.ip(1)
                        warning('The ip of %d-th well is not valid, please check its data range.', k);
                    end
                case 'rho' 
                    if maxVal > options.range.rho(2) || minVal < options.range.rho(1)
                        warning('The rho of %d-th well is not valid, please check its data range.', k);
                    end
                case 'vpvs_ratio' 
                    if maxVal > options.range.vpvs_ratio(2) || minVal < options.range.vpvs_ratio(1)
                        warning('The vpvs_ratio of %d-th well is not valid, please check its data range.', k);
                    end
                case 'possion' 
                    if maxVal > options.range.possion(2) || minVal < options.range.possion(1)
                        warning('The possion of %d-th well is not valid, please check its data range.', k);
                    end
            end
        end
        
        
    end
    
end