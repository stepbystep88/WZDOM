function [index, returnHeader] = bsIndexOfTraceSetOnInIdAndCrossId(GSegyInfo, inId, crossId, index)
% find the trace index of a segy file by inline id and crossline id
%
% Programmed by: Bin She (Email: bin.stepbystep@gmail.com)
% Programming dates: Nov 2019
% -------------------------------------------------------------------------
% Input
% GSegyInfo.fid       file handle
% GSegyInfo.volHeader volume header of the segy file
% inId      inline id
% crossId   crossline id
%
% Output
% index     the index of a trace matching the input inline and crossline ids
    sizeTrace = GSegyInfo.volHeader.sizeTrace+240;
    returnHeader = [];
    
    % check whether index work
    if nargin > 3 && index <= GSegyInfo.volHeader.traceNum && index > 0
        offset = 3600 + sizeTrace*(index-1);    
        fseek(GSegyInfo.fid, offset, -1);                 
        returnHeader = bsReadTraceHeader(GSegyInfo);
        if(returnHeader.crossId==crossId && returnHeader.inId==inId)
            return;
        else
            returnHeader = [];
        end

    end
    
    fseek(GSegyInfo.fid, 3600, -1);       % skip 3600 bytes
    
    startPos = 1;
    endPos = GSegyInfo.volHeader.traceNum;
    index = floor((startPos+endPos) / 2);

    while(startPos~=index && index~=endPos)
        % bisearch
        index = floor((startPos+endPos) / 2);

        offset = 3600 + sizeTrace*(index-1);    
        fseek(GSegyInfo.fid, offset, -1);                 
        trHeader = bsReadTraceHeader(GSegyInfo);    

        if( trHeader.inId > inId)
            endPos = index - 1;
        elseif( trHeader.inId < inId)
            startPos = index + 1;
        else 
            if(trHeader.crossId > crossId)
                endPos = index - 1;
            elseif(trHeader.crossId < crossId)
                startPos = index + 1;
            else
                break;
            end
        end
    end

    if(trHeader.inId ~= inId || trHeader.crossId ~= crossId)
        % no such trace matching the given inline and crossline ids
        index = -1;
    else
        % seek to the current index
        fseek(GSegyInfo.fid, 3600+sizeTrace*(index-1), -1);

        while(true)
            trHeader = bsReadTraceHeader(GSegyInfo);    % read trace header
            if(trHeader.crossId==crossId && trHeader.inId==inId)
                index = index - 1;
                if( index == 0)
                    break;
                else
                    % seek to the previous trace until it doesn't match the
                    % input inline and crossline ids
                    % this while loop ensures the output index is the first
                    % trace satisfying the criterion
                    fseek(GSegyInfo.fid, -sizeTrace-240, 0);
                end
            else
                break;
            end
        end
        
        index = index + 1;
    end
end