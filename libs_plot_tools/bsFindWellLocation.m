function [wellPos, wellIndex, wellNames] = bsFindWellLocation(wellLogs, inIds, crossIds)

    wellPos = [];
    wellIndex = [];
    wellNames = {};
    
    if isempty(wellLogs)
        return;
    end
    
    wells = cell2mat(wellLogs);
    wellInIds = [wells.inline];
    wellCrossIds = [wells.crossline];
    
    for i = 1 : length(inIds)
        for j = 1 : length(wellInIds)
            if wellInIds(j) == inIds(i) && wellCrossIds(j) == crossIds(i)
                wellPos = [wellPos, i];
                wellIndex = [wellIndex, j];
                wellNames = [wellNames, wellLogs{j}.name];
            end
        end
    end
end