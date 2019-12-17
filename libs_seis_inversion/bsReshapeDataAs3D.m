function volume = bsReshapeDataAs3D(data, nInline, nCrossline)
    [sampNum, nTrace] = size(data);
    assert(nTrace==nInline*nCrossline, 'the number of traces must be equal to #inline * #crossline');
    
    volume = reshape(data, sampNum, nCrossline, nInline);
    
    % let the volume follows the order of sampNum, inline, crossline
    volume = permute(volume, [1 3 2]);

    
end