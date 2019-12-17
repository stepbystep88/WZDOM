function [volume, rangeInline, rangeCrossline, startTime] ...
    = bsReadTracesAs3DVolume(fileName, GSegyInfo, varargin)
%% read volume from a segy file
%
% Programmed by: Bin She (Email: bin.stepbystep@gmail.com)
% Programming dates: Dec 2019
% -------------------------------------------------------------------------

    p = inputParser;
    
    [rangeInline, rangeCrossline] = bsGetWorkAreaRange(GSegyInfo, fileName);
    
    addParameter(p, 'rangeInline', rangeInline);
    addParameter(p, 'rangeCrossline', rangeCrossline);
    addParameter(p, 'timeLine', []);
    addParameter(p, 'dt', []);
    addParameter(p, 'upNum', []);
    addParameter(p, 'downNum', []);
    addParameter(p, 't0', []);
    
    p.parse(varargin{:});  
    options = p.Results;
    
    rangeInline = options.rangeInline;
    rangeCrossline = options.rangeCrossline;
    
    [inIds, crossIds] = bsGetCDPsByRange(rangeInline, rangeCrossline);
    
    if isempty(options.timeLine)
        
        if ~isempty(options.t0)
            pos = floor((options.t0 - GSegyInfo.t0)/options.dt);
            t0 = pos*options.dt + GSegyInfo.t0;
            startTime = t0 * ones(1, length(inIds));
            
        else
            startTime = GSegyInfo.t0 * ones(1, length(inIds));
            pos = 0;
        end
        data = bsReadTracesByIds(fileName, GSegyInfo, inIds, crossIds);
        data = data(pos+1:end, :);
        
    else
        dt = options.dt;
        horizonTime = bsGetHorizonTime(options.timeLine, inIds, crossIds, 1);
        startTime = horizonTime - options.upNum*dt;
        sampNum = options.upNum + options.downNum;
        
        data = bsReadTracesByIds(fileName, GSegyInfo, inIds, crossIds, startTime, sampNum, dt);
    end
    
    nInline = rangeInline(end) - rangeInline(1) + 1;
    nCrossline = rangeCrossline(end) - rangeCrossline(1) + 1;
    
    volume = bsReshapeDataAs3D(data, nInline, nCrossline);
    startTime = reshape(startTime, nCrossline, nInline);
    startTime = permute(startTime, [2, 1]);
end