function algParams = bsCalNewValueOfAlgParams(algParams)
    
    
    names = fieldnames(algParams);
    nField = length(names);
    
    for i = 1 : nField
        t = getfield(algParams, names{i});
        
        nTrain = t.nTrain;
        
%         if nTrain <= 1
%             continue;
%         end
          
        if nTrain >= size(t.store, 1)
            

            sortData = sortrows(t.store(1:nTrain, :), -2);
            value = mean(sortData(1:5, 1));
%             sortData = sortData(1:10, :);
%             weight = sortData(:, 2) ./ sum(sortData(:, 2));
%             S = sortData(:, 1);
%             value = sum(weight .* S .^ 2) / sum(weight .* S);
%             
%             trueVal = normrnd(value, 0.1);
            trueVal = value + randi([-0.1,0.1]*10000)/10000;
%             trueVal = value;

            if trueVal > t.max
                t.value = t.max;
            elseif trueVal < t.min
                t.value = t.min;
            else
                t.value = trueVal;
            end
            
%             if rand() < 0.1
%                 t.value = t.min + randi([1, floor(10000*(t.max-t.min))])/10000;
%             end
        else
            t.value = t.min + randi([1, floor(10000*(t.max-t.min))])/10000;
        end
        
        
        algParams = setfield(algParams, names{i}, t);
    end
end