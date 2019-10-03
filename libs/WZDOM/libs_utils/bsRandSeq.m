function [ out ] = bsRandSeq(data, k)
%% randomly generate a sequence of length k from data

    n = length(data);
    
    for i = 1 : k
        index = floor( (n-i+1)*rand() ) + i;
       
        
        temp = data(index);
        data(index) = data(i);
        data(i) = temp;
    end

    out = data(1:k);
end

