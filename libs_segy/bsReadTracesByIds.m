function [trData, GSegyInfo, trHeader] = bsReadTracesByIds(fileName, GSegyInfo, inIds, crossIds, startPos, sampNum)
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
% startPos      the start position of a trace where we extract the data from
% sampNum       how many sample points we would extract from a trace
% -------------------------------------------------------------------------
            
    GSegyInfo = bsReadVolHeader(fileName, GSegyInfo);        

    [~, trNum] = size(inIds);                            
    
    if nargin > 4
        trData = zeros(sampNum, trNum);
    else
        trData = zeros(GSegyInfo.volHeader.sampNum, trNum);
    end

%     volHeader.dataForm = 5;
    
    trHeader = [];
    for i = 1 : trNum
        
        inId = inIds(i);
        crossId = crossIds(i);
        
        index = bsIndexOfTraceSetOnInIdAndCrossId(GSegyInfo, inId, crossId);
        fseek(GSegyInfo.fid, 3600 + (index-1)*(240+GSegyInfo.volHeader.sizeTrace), -1);
        
        trHeader = bsReadTraceHeader(GSegyInfo);
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
            trData(:, i) = data(sPos : ePos);
        end

        
    end

    fclose(GSegyInfo.fid);                                        
end