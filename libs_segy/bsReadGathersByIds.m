function gathers = bsReadGathersByIds(fileName, GSegyInfo, inIds, crossIds, startTime, sampNum, dt)
%% read gathers from a segy file with given inline and crossline ids
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
    if nargin > 4
        startPos = bsCalcT0Pos(GSegyInfo, startTime, dt);
    end
    
    [~, pointNum] = size(inIds);
    gathers = cell(1, pointNum);
    
    for i = 1 : pointNum
        inId = inIds(i);
        crossId = crossIds(i);
        index = bsIndexOfTraceSetOnInIdAndCrossId(GSegyInfo, inId, crossId);
        
        data = [];
        offsets = [];
        
        if(index == -1)
            warning("The trace of inline=%d and cossline=%d doesn't exist in the file %s!\n", ...
                inId, crossId, fileName);
            continue;
        end
        
        % skip to the first trace matching inId and crossId
        fseek(GSegyInfo.fid, 3600 + (index-1)*(240+GSegyInfo.volHeader.sizeTrace), -1);
        
        while (index <= GSegyInfo.volHeader.traceNum)
            index = index + 1;
            trHeader = bsReadTraceHeader(GSegyInfo);
            trData = bsReadTraceData(GSegyInfo);
            
            % this is the next gather
            if (trHeader.inId ~= inId || trHeader.crossId ~= crossId)  
                break;
            end
            
            offsets = [offsets, trHeader.offset];
            if nargin > 4
                iPos = startPos(i);
                trData = trData(iPos+1 : iPos+sampNum);
                data = [data, trData];
            else
                data = [data, trData];
            end
            
        end
        
        gathers{i}.data = data;
        gathers{i}.offsets = offsets;
    end
    
    fclose(GSegyInfo.fid);                                        % 读取完毕之后需要关闭fin
    
end