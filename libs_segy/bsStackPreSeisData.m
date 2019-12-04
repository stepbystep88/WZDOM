function postSeisData = bsStackPreSeisData(fileName, GSegyInfo, inIds, crossIds, startTime, sampNum, dt)
%% read and calculate post-stack seismic data from a pre-stack segy file with given inline and crossline ids
%
% Programmed by: Bin She (Email: bin.stepbystep@gmail.com)
% Programming dates: Nov 2019
% -------------------------------------------------------------------------
    trNum = length(inIds);
    postSeisData = zeros(sampNum, trNum);
    gathers = bsReadGathersByIds(fileName, GSegyInfo, inIds, crossIds, startTime, sampNum, dt);
    
    for i = 1 : trNum
        if ~isempty(gathers{i}.data)
            postSeisData(:, i) = mean(gathers{i}.data, 2);
        end
    end
    
end
