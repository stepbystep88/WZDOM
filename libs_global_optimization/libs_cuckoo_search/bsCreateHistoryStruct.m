function [parameters] = bsCreateHistoryStruct(paramsConfigure, nHistory)
    nParameter = size(paramsConfigure, 1);
    
    parameters = [];
    
    for i = 1 : nParameter
        
        t.index = 0;
        t.nTrain = 0;
        t.store = zeros(nHistory, 2);
        
        t.min = paramsConfigure{i, 2};
        t.max = paramsConfigure{i, 3};
        
%         t.store(1, 1) = paramsConfigure{i, 1};
        
        t.value = paramsConfigure{i, 1};
        
        parameters = setfield(parameters, paramsConfigure{i, 4}, t);
    end
    
end