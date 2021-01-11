function wellLogs = bsResetValidWellData(GInvParam, wellLogs, varargin)
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
            data = wellLogs{k}.wellLog(:, index);
            maxVal = max(data);
            minVal = min(data);
            
            if isempty(maxVal) || isempty(minVal)
                warning('The %s welllog data of %d-th well is empty', name_i, k);
                continue;
            end
            invalid_index = [];
            switch  name_i
                case 'vp' 
                    invalid_index = find(data > options.range.vp(2) | data < options.range.vp(1));
                case 'vs' 
                    invalid_index = find(data > options.range.vs(2) | data < options.range.vs(1));
                case 'ip' 
                    invalid_index = find(data > options.range.ip(2) | data < options.range.ip(1));
                case 'rho' 
                    invalid_index = find(data > options.range.rho(2) | data < options.range.rho(1));
                case 'vpvs_ratio' 
                    invalid_index = find(data > options.range.vpvs_ratio(2) | data < options.range.vpvs_ratio(1));
                case 'possion' 
                    invalid_index = find(data > options.range.possion(2) | data < options.range.possion(1));
            end
            
            if  ~isempty(invalid_index)
                for j = 1 : length(invalid_index)
                    if invalid_index(j) == 1
                        wellLogs{k}.wellLog(invalid_index(j), index) = 0.5*(wellLogs{k}.wellLog(2, index) + wellLogs{k}.wellLog(3, index));
                    elseif invalid_index(j)==length(data)
                        wellLogs{k}.wellLog(invalid_index(j), index) = 0.5*(wellLogs{k}.wellLog(length(data)-1, index) + wellLogs{k}.wellLog(length(data)-2, index));
                    else
                        wellLogs{k}.wellLog(invalid_index(j), index) = 0.5*(wellLogs{k}.wellLog(invalid_index(j)-1, index) + wellLogs{k}.wellLog(invalid_index(j)+1, index));
                    end
                    
                end
            end
        end
        
        
    end
    
end