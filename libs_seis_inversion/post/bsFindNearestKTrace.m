function ids = bsFindNearestKTrace(iTrace, xs, ys, K, nTracePerLine)
    sPos = iTrace - K - nTracePerLine * K;
    if (sPos < 1)
        sPos = 1;
    end
    
    ePos = iTrace + K + nTracePerLine * K;
    if ePos > length(xs)
        ePos = length(xs);
    end
    
    subIndex = sPos : ePos;
    dist = (xs(iTrace) - xs(subIndex)).^2 + (ys(iTrace) - ys(subIndex)).^2;
    
    [~, index] = bsMinK([dist', subIndex'], K);
    ids = subIndex(index);
end