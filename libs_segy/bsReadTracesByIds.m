function [trData, GSegyInfo] = bsReadTracesByIds(fileName, GSegyInfo, inIds, crossIds, startTime, sampNum, dt)
%% read traces from a segy file with given inline and crossline ids
%
% Programmed by: Bin She (Email: bin.stepbystep@gmail.com)
% Programming dates: Nov 2019
% -------------------------------------------------------------------------
% Input
% fileName      file name of the input segy file
% GSegyInfo     basical infomation of the segy file
% inIds         inline ids
% crossIds      crossline ids
% startTime     the start time of a trace that we extract the data from
% sampNum       how many sample points we would extract from a trace
% -------------------------------------------------------------------------
            
    GSegyInfo = bsReadVolHeader(fileName, GSegyInfo);        
    sizeTrace = GSegyInfo.volHeader.sizeTrace+240;    

    [~, trNum] = size(inIds);                            
    
    if nargin > 4
        trData = zeros(sampNum, trNum);
        % the start position of a trace that we extract the data from
        startPos = bsGetT0Pos(GSegyInfo, startTime, dt);
    else
        trData = zeros(GSegyInfo.volHeader.sampNum, trNum);
    end

%     volHeader.dataForm = 5;
    
    for i = 1 : trNum
        if mod(i, 10000) == 0
            % print information
            fprintf('Reading %d%% data from segy file %s...\n', round(i/trNum*100), fileName);
        end
        
        if i == 1
            index = bsIndexOfTraceSetOnInIdAndCrossId(GSegyInfo, inIds(i), crossIds(i));
        else
            index = bsIndexOfTraceSetOnInIdAndCrossId(GSegyInfo, inIds(i), crossIds(i), index);
        end
        
        if index > 0
            fseek(GSegyInfo.fid, 3600+sizeTrace*(index-1), -1);
            trHeader = bsReadTraceHeader(GSegyInfo);
            index = index + 1;
        else
            warning('Trace inline=%d, crossline=%d can not be found in file %s', inIds(i), crossIds(i), fileName);
            continue;
        end
        
        data = bsReadTraceData(GSegyInfo);

        if nargin > 4
            iPos = startPos(i);
            
            % just in case the original data is too short to extract the
            % target part
            if iPos < 0
                sPos = 1;
                lsPos = abs(iPos) + 1;
            else
                sPos = iPos + 1;
                lsPos = 1;
            end
            
            if iPos + sampNum > length(data)
                ePos = length(data);
                lePos = length(data) - iPos;
            else
                ePos = iPos+sampNum;
                lePos = sampNum;
            end
            trData(lsPos : lePos, i) = data(sPos : ePos);
        else
            trData(:, i) = data;
        end

        
    end

    fclose(GSegyInfo.fid);    
    

end