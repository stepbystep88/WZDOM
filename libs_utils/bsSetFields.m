function data = bsSetFields(data, pairs)
%% set multiple fields of a struct
% Bin She, bin.stepbystep@gmail.com, Nov, 2019
% Programmed by: Bin She (Email: bin.stepbystep@gmail.com)
% Programming dates: May 2019
% -------------------------------------------------------------------------

    
    nFields = size(pairs, 1);
    
    for i = 1 : nFields
        data = setfield(data, pairs{i, 1}, pairs{i, 2});
    end
    
    
    
end