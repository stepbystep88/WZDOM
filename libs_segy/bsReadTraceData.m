function [traceData] = bsReadTraceData(GSegyInfo)
%% write volume header into a segy
%
% Programmed by: Bin She (Email: bin.stepbystep@gmail.com)
% Programming dates: Nov 2019
% -------------------------------------------------------------------------
% Input
% GSegyInfo         a struct containing some basical information of a
% related segy file
%
% Output
% traceData         the data of current trace
% -------------------------------------------------------------------------
    
    switch GSegyInfo.volHeader.dataForm
        case 5 % IEEE format
            traceData = fread(GSegyInfo.fid, GSegyInfo.volHeader.sampNum, 'float32');             
        otherwise % read as IBM format
            traceData = fread(GSegyInfo.fid, GSegyInfo.volHeader.sampNum, 'uint32=>uint32');      
            traceData = ibm2double(traceData);                      % IBM_To_Double
    end

    if GSegyInfo.isNegative
        traceData = -traceData;
    end
end