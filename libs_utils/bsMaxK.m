function [B, index] = bsMaxK(A, k, dim)
    
    if size(A, 1) == 1
        A = A';
    end
    
    if ~exist('dim', 'var')
        dim = -1;
    end
        
    [B, I] = sortrows(A, dim);
    B = B(1:k, :);
    index = I(1:k);
end