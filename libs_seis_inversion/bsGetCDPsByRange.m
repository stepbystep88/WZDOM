function [inIds, crossIds] = bsGetCDPsByRange(rangeInline, rangeCrossline)

    [X, Y] = meshgrid(rangeInline(1):rangeInline(end), rangeCrossline(1):rangeCrossline(end));
    [n1, n2] = size(X);
    
    inIds = reshape(X, 1, n1 * n2);
    crossIds = reshape(Y, 1, n1 * n2);
end