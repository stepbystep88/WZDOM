function timeLines = bsRestoreHorizon(timeLines)
    for k = 1 : length(timeLines)
        timeLine = timeLines{k};
        
        % sortTimeLine
        timeLine = sortrows(timeLine, [1, 2]);
        
        firstCrossline = min(unique(timeLine(:, 2)));
        endCrossline = max(unique(timeLine(:, 2)));
        firstInline = min(unique(timeLine(:, 1)));
        endInline = max(unique(timeLine(:, 1)));
        
        nInline = endInline - firstInline + 1;
        nCrossline = endCrossline - firstCrossline + 1;
        nTrace = nInline * nCrossline;
        
        [X, Y] = meshgrid(firstInline : endInline, firstCrossline : endCrossline);
        Z = zeros(size(X));
        
        parfor i = 1 : nInline
            for j = 1 : nCrossline
                inId = X(j, i);
                crossId = Y(j, i);
                
               
                [~, ~, Z(j, i)] = bsGetWellBaseInfo(timeLine, ...
                    inId, crossId, 1, 2, 1, 2, 3);
                
                
            end
            if mod(i, 5) == 0
                fprintf('Calculating horizon information %d/%d\n', i, nInline);
            end
        end
        
        newTimeLine = [reshape(X, nTrace, 1), reshape(Y, nTrace, 1), reshape(Z, nTrace, 1)];
        timeLines{k} = newTimeLine;
        
    end
end