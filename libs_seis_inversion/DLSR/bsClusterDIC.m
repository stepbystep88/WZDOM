function [newDIC, newIdx] = bsClusterDIC(DIC, K) 
    Idx =kmeans(DIC', K, 'Replicates', 10);
    newDIC = zeros(size(DIC));
    j = 1;
    
    newIdx = zeros(1, size(DIC, 2));
    
    for i = 1 : K
        index = find(Idx == i);
        newDIC(:, j:j+length(index)-1) = DIC(:, index);
        newIdx(j:j+length(index)-1) = i;
        j = j + length(index);
        
        
    end
end