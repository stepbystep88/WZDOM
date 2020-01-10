function wellLogs = bsSetNameForWelllogs(wellLogs)
    for i = 1 : length(wellLogs)
        wellLogs{i}.name = sprintf('W%d',i);
    end
end