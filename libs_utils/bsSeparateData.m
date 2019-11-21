function [cdata, idxs] = bsSeparateData(data, nGroup)
    nCase = length(data);
    for i = 1 : nCase-1
        assert(length(data{i}) == length(data{i+1}), 'The length of input data must be the same!');
    end
    
    nData = length(data{1});
    nDataInGroup = floor(nData/nGroup);
    cdata = cell(nCase, nGroup);
    idxs = cell(1, nGroup);
    
    for i = 1 : nGroup
        spos = nDataInGroup * (i - 1) + 1;
        
        if i == nGroup
            epos = nData;
        else
            epos = nDataInGroup * i;
        end
        
        for j = 1 : nCase
            cdata{j, i} = data{j}(spos:epos);
        end
        
        idxs{i} = spos : epos;
    end
end