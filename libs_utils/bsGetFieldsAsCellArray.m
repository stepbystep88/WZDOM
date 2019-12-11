function [values] = bsGetFieldsAsCellArray(data)
%% get the value of input field names in pairs
% Bin She, bin.stepbystep@gmail.com, Dec, 2019
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    fields = fieldnames(data);
    nField = length(fields);
    
    values = cell(1, 2*nField);
    
    for i = 1 : nField
        values{i*2-1} = fields{i};
        values{i*2} = getfield(data, fields{i});
    end
    
end