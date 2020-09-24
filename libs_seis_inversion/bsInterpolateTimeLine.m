function newLine = bsInterpolateTimeLine(timeLine)

    nLine = length(timeLine);
    rangeInline = [min(timeLine{1}(:, 1)), max(timeLine{1}(:, 1))];
    rangeCrossline = [min(timeLine{1}(:, 2)), max(timeLine{1}(:, 2))];
    
    nInline = rangeInline(2) - rangeInline(1) + 1;
    nCrossline = rangeCrossline(2) - rangeCrossline(1) + 1;
    
    newLine = cell(1, nLine);
    for i = 1 : nLine
        iLine = timeLine{i};
        if nInline * nCrossline == size(iLine, 1)
            newLine{i} = iLine;
            continue;
        end
        
        [xq, yq] = meshgrid(rangeInline(1):rangeInline(2), rangeCrossline(1):rangeCrossline(2));
        
        x = iLine(:, 1);
        y = iLine(:, 2);
        v = iLine(:, 3);
        F = scatteredInterpolant(x, y, v);
        F.Method = 'natural';
        vq = F(xq, yq);
        
        newLine{i} = [reshape(xq, [], 1), reshape(yq, [], 1), reshape(vq, [], 1)];
    end
    
end