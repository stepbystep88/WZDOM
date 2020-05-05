function ids = bsFindNearestKTrace(iTrace, xs, ys, K, nTracePerLine)
    sPos = iTrace - K - nTracePerLine * K;
    if (sPos < 1)
        sPos = 1;
    end
    
    ePos = iTrace + K + nTracePerLine * K;
    if ePos > length(xs)
        ePos = length(xs);
    end
    
    dist = (xs(iTrace) - xs(sPos:ePos)).^2 + (ys(iTrace) - ys(sPos:ePos)).^2;
    
    [~, ids] = bsMinK(dist, K);
end