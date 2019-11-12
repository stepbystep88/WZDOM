function [trData, inIds, crossIds, trHeaders] = bsReadAllTraces(fileName, GSegyInfo, isSaveHeaders)
%% write volume header into a segy
%
% Programmed by: Bin She (Email: bin.stepbystep@gmail.com)
% Programming dates: Nov 2019
% -------------------------------------------------------------------------
% Input
% fileName          file name of a segy file
% GSegyInfo         a struct containing some basical information of the
% segy file, it will also be output 
% isSaveHeaders     whether return the trace headers of the segy file,
% yes -> trHeaders will be a cell, no -> trHeaders will be the trace header
% of the last trace
%
% Output
% volHeader         volume header
% GSegyInfo         a struct containing some basical information of the segy file
% trHeaders         see isSaveHeaders
% inIds         
% crossIds

% -------------------------------------------------------------------------

    if nargin == 2
        isSaveHeaders = 0;
    end

    GSegyInfo = bsReadVolHeader(fileName, GSegyInfo);        
%     volHeader.sampNum = 180;
    
    volHeader = GSegyInfo.volHeader;
    
    trNum = volHeader.traceNum;                  
    trData = zeros(volHeader.sampNum, trNum);
    inIds = zeros(1, trNum);
    crossIds = zeros(1, trNum);
    
    if isSaveHeaders
        trHeaders = cell(1, trNum);
    end
            
    % read data trace by trace
    for i = 1 : trNum
        
        trHeader = bsReadTraceHeader(GSegyInfo);
        
        if isSaveHeaders
            trHeaders{i} = trHeader;
        else
            trHeaders = trHeader;
        end
        
        
        inIds(i) = trHeader.inId;
        crossIds(i) = trHeader.crossId;
        
%         fseek(GSegyInfo.fid, 240, 0);
        data = bsReadTraceData(GSegyInfo);
        trData(:, i) = data;   
    end
  
    
    fclose(GSegyInfo.fid);                 
end