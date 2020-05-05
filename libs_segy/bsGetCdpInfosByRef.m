function [infos] = bsGetCdpInfosByRef(fileName, GSegyInfo)
    GSegyInfo = bsReadVolHeader(fileName, GSegyInfo);        
%     volHeader.sampNum = 180;
    
    volHeader = GSegyInfo.volHeader;
    
    trNum = volHeader.traceNum;            

    infos = zeros(trNum, 4);
    
    % read data trace by trace
    for i = 1 : trNum
        trHeader = bsReadTraceHeader(GSegyInfo);
       
        infos(i, :) = [trHeader.inId, trHeader.crossId, trHeader.X, trHeader.Y];
        
        fseek(GSegyInfo.fid, volHeader.sizeTrace, 0);
    end
    
    fclose(GSegyInfo.fid);
end