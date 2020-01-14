function volume = bsReshapeDataAs2D(data)
    % the original data follows the order of sampNum, inline, crossline, so
    % we need to premute the data
    data = permute(data, [1 3 2]);
    
    [sampNum, nCrossline, nInline] = size(data);
    volume = reshape(data, sampNum, nCrossline*nInline);
end