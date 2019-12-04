function [inIds, crossIds] = bsGetProfileCrossingWells(GInvParam, wellLogs, varargin)

    p = inputParser;
    
    [rangeInline, rangeCrossline] = bsGetWorkAreaRangeByParam(GInvParam);
    
    addParameter(p, 'isAlongCrossline', 1);
    addParameter(p, 'method', 'cubic');
    addParameter(p, 'rangeInline', rangeInline);
    addParameter(p, 'rangeCrossline', rangeCrossline);
    
    p.parse(varargin{:});  
    options = p.Results;
    
    wells = cell2mat(wellLogs);
    wellInIds = [wells.inline];
    wellCrossIds = [wells.crossline];
    nWell = length(wellLogs);
    
    % get the range of current work area
    rangeInline = options.rangeInline;
    rangeCrossline = options.rangeCrossline;
    
    if nWell == 1
        if options.isAlongCrossline
            traceNum = rangeCrossline(2) - rangeCrossline(1) + 1;
            
            inIds = ones(1, traceNum) * wellInIds(1);
            crossIds = rangeCrossline(1) : rangeCrossline(2);
        else
            
            traceNum = rangeInline(2) - rangeInline(1) + 1;
            
            crossIds = ones(1, traceNum) * wellCrossIds(1);
            inIds = rangeInline(1) : rangeInline(2);
        end
    else
        if options.isAlongCrossline
            [inIds, crossIds] = bsInterpolateALine(wellInIds, wellCrossIds, ...
                rangeInline, rangeCrossline, options.method);
        else
            [crossIds, inIds] = bsInterpolateALine(wellCrossIds, wellInIds, ...
                rangeCrossline, rangeInline, options.method);
        end
    end
end

function [outInIds, outCrossIds] = bsInterpolateALine(inIds, crossIds, ...
    rangeInline, rangeCrossline, interp_method)

    setInIds = [inIds(1), inIds, inIds(length(inIds))];
    setCrossIds = [rangeCrossline(1), crossIds, rangeCrossline(2)];
    
    
    outCrossIds = rangeCrossline(1) : 1 : rangeCrossline(2);
    outInIds = interp1(setCrossIds, setInIds, outCrossIds, interp_method);

    for i = 1 : length(outInIds)
        outInIds(i) = floor(outInIds(i));
    end
    
    for i = 1 : length(crossIds)
        index = outCrossIds == crossIds(i);
        outInIds(index) = inIds(i);
    end
    
    outInIds(outInIds > rangeInline(2)) = rangeInline(2);
    outInIds(outInIds < rangeInline(1)) = rangeInline(1);
end