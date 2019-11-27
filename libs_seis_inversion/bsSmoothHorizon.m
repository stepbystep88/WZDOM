function timeLines = bsSmoothHorizon(timeLines, N)
    for i = 1 : length(timeLines)
        timeLine = timeLines{i};
        [inRange, crossRange] = bsGetRange(timeLine);
        
        [timeSlice, inGrid, crossGrid] = bsLine2Slice(timeLine, inRange, crossRange);
        
%         h = fspecial('gaussian', [3, 3]);
        smoothTimeSlice = smooth2a(timeSlice, N, N);
        
        nData = length(smoothTimeSlice(:));
        
        smoothTimeLine = [reshape(inGrid, nData, 1), ...
            reshape(crossGrid, nData, 1), ...
            reshape(smoothTimeSlice, nData, 1)];
        
        timeLines{i} = smoothTimeLine;
    end
end

function [inRange, crossRange] = bsGetRange(timeLine)
    inRange = [min(timeLine(:, 1)), max(timeLine(:, 1))];
    crossRange = [min(timeLine(:, 2)), max(timeLine(:, 2))];
end

function [timeSlice, inGrid, crossGrid] = bsLine2Slice(timeLine, inRange, crossRange)
    [inRange, crossRange] = bsGetRange(timeLine);
    [inGrid, crossGrid] = meshgrid(inRange(1):inRange(2), crossRange(1):crossRange(2));
    
    timeSlice = zeros(size(inGrid));
    
    for i = 1 : size(timeLine, 1)
        inId = timeLine(i, 1);
        crossId = timeLine(i, 2);
        
        y = inId - inRange(1) + 1;
        x = crossId - crossRange(1) + 1;
        timeSlice(x, y) = timeLine(i, 3);
    end
    
    [x, y] = find(timeSlice == 0);
    [xx, yy] = find(timeSlice ~= 0);
    
    for i = 1 : length(x)
        
        dist = (xx-x(i)).^2 + (yy-y(i)).^2;
        [~, index] = min(dist);
        
        timeSlice(x(i), y(i)) = timeSlice(xx(index), yy(index));
    end
end
