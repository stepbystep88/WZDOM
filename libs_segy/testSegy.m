%% test read and write funciton of this segy toolbox
inFileName = 'in.sgy';
outFileName = 'out.sgy';

load GSegyInfo.mat;
GInSegyInfo = bsReadVolHeader(inFileName, GSegyInfo);
GOutSegyInfo = bsWriteVolHeader(outFileName, GInSegyInfo);

inData = zeros(GInSegyInfo.volHeader.sampNum, GInSegyInfo.volHeader.traceNum);

for i = 1 : GInSegyInfo.volHeader.traceNum
    trHeader = bsReadTraceHeader(GInSegyInfo);
    data = bsReadTraceData(GInSegyInfo);
    inData(:, i) = data;
    
    bsWriteTrace(GOutSegyInfo, trHeader, -data);
end

fclose(GInSegyInfo.fid);
fclose(GOutSegyInfo.fid);

[outData, inIds, crossIds, trHeaders] = bsReadAllTraces(outFileName, GSegyInfo);

figure;
subplot(2, 1, 1);   imagesc(inData);
subplot(2, 1, 2);   imagesc(outData);