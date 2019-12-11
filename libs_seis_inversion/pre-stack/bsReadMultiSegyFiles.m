function seisData = bsReadMultiSegyFiles(separates, inline, crossline, startTime, sampNum, dt)
    nFile = length(separates);
    seisData = zeros(sampNum, nFile);
    for i = 1 : nFile
        separate = separates(i);
        seisData(:, i) = bsReadTracesByIds(separate.fileName, separate.segyInfo, inline, crossline, startTime, sampNum, dt);
    end
end