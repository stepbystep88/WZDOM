function [wellLogs, timeLine] = bsPostGenWelllogs(GInvParam, trueModel, well_ids, t0)
    nWell = length(well_ids);
    
    wellLogs = cell(1, nWell);
    [sampNum, traceNum] = size(trueModel);
    
    time = (0:sampNum-1)'*GInvParam.dt + t0;
    for i = 1 : nWell
        t.name = sprintf('#%d', i);
        t.inline = 1;
        t.crossline = well_ids(:, i);
        t.wellLog = [trueModel(:, well_ids(:, i)), time];
        
        wellLogs{i} = t;
    end
    
    timeLine = {ones(1, traceNum) * t0};
end