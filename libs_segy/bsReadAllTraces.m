function [trData, inIds, crossIds] = bsReadAllTraces(fileName, GSegyInfo, startTime, sampNum, dt)
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


    GSegyInfo = bsReadVolHeader(fileName, GSegyInfo);        
%     volHeader.sampNum = 180;
    
    volHeader = GSegyInfo.volHeader;
    
    trNum = volHeader.traceNum;            
    if nargin == 5
        trData = zeros(sampNum, trNum);
        % the start position of a trace that we extract the data from
        startPos = bsCalcT0Pos(GSegyInfo, startTime, dt);
    else
        trData = zeros(volHeader.sampNum, trNum);
    end
    
    
    inIds = zeros(1, trNum);
    crossIds = zeros(1, trNum);
    
    % read data trace by trace
    for i = 1 : trNum
        if mod(i, 10000) == 0
            % print information
            fprintf('Reading %d%% data from segy file %s...\n', round(i/trNum*100), fileName);
        end
        
        trHeader = bsReadTraceHeader(GSegyInfo);
 
        inIds(i) = trHeader.inId;
        crossIds(i) = trHeader.crossId;
        
%         fseek(GSegyInfo.fid, 240, 0);
        data = bsReadTraceData(GSegyInfo);
        
        if nargin == 5
            iPos = startPos(i);
            trData(:, i) = data(iPos+1 : iPos+sampNum);
        else
            trData(:, i) = data;
        end
    end
  
    
    fclose(GSegyInfo.fid);                 
end