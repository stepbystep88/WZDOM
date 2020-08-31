function data = bsInterpolate3DData(nTrace, wellData, weights, indexies)

    sampNum = size(wellData, 1);
    data = zeros(sampNum, nTrace);
    
    for i = 1 : nTrace
        weight = weights(:, i);
        usedWellData = wellData(:, indexies(i, :));
        
        V = usedWellData * weight;
        
        data(:, i) = V;
        
        if mod(i, 10000) == 0
            fprintf('Interpolating the model of trace %d/%d...\n', i, nTrace);
        end
    end
end