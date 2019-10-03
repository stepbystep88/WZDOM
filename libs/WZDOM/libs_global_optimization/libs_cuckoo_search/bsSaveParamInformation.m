function history = bsSaveParamInformation(history, algParams, iter, maxIter)
    
    
    names = fieldnames(algParams);
    nField = length(names);
    
    if iter == 1
        history = zeros(nField, maxIter);
    end
    
    for i = 1 : nField
        t = getfield(algParams, names{i});
        history(i, iter) = t.value;
    end
end