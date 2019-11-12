function [index] = bsIndexOf(array, value)
% find the index of a value in an array
%
% Programmed by: Bin She (Email: bin.stepbystep@gmail.com)
% Programming dates: Nov 2019
% -------------------------------------------------------------------------
% Input
% array     an array
% value     the entry to be found
%
% Output
% index     the index of the entry 

    [~, col] = size(array);       

    startPos = 1;
    endPos = col;
    index = floor((startPos+endPos) / 2);
    maxError = 0.00001;

    while(startPos~=index && index~=endPos)
        index = floor((startPos+endPos) / 2);

        if( array(index) - value > maxError)
            endPos = index-1;
        elseif( value - array(index) > maxError)
            startPos = index+1;
        else
            break;
        end
    end

end